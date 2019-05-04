#include "multiplane.h"

#include "util/error.h"
#include <algorithm>
#include <iostream>
#include <thrust/device_vector.h>

Multiplane MultiplaneBuilder::getCuMultiPlane() {
    std::vector<CompositeLens> data;
    std::vector<int> subl_size;

    for (size_t i = 0; i < m_builders.size(); i++) {
        subl_size.push_back(m_builders[i].length());
        auto lens = m_builders[i].getCuLens();
        data.push_back(lens);
    }

    if (data.size() == 0 || m_source_z.size() == 0) {
        std::cerr << "No lens and/or source planes given " << data.size() << "-"
                  << m_source_z.size() << std::endl;
        throw(-1);
    }

    prepare();

    CompositeLens *lens_ptr;
    float *src_ptr, *dist_lens_ptr, *dist_src_ptr;

    size_t lens_size = sizeof(CompositeLens) * data.size();
    gpuErrchk(cudaMalloc(&lens_ptr, lens_size));
    gpuErrchk(
        cudaMemcpy(lens_ptr, &data[0], lens_size, cudaMemcpyHostToDevice));

    size_t src_size = sizeof(float) * m_source_z.size();
    gpuErrchk(cudaMalloc(&src_ptr, src_size));
    gpuErrchk(
        cudaMemcpy(src_ptr, &m_source_z[0], src_size, cudaMemcpyHostToDevice));

    size_t dist_lens_size = sizeof(float) * m_dists_lenses.size();
    gpuErrchk(cudaMalloc(&dist_lens_ptr, dist_lens_size));
    gpuErrchk(cudaMemcpy(dist_lens_ptr, &m_dists_lenses[0], dist_lens_size,
                         cudaMemcpyHostToDevice));

    size_t dist_src_size = sizeof(float) * m_dists_sources.size();
    gpuErrchk(cudaMalloc(&dist_src_ptr, dist_src_size));
    gpuErrchk(cudaMemcpy(dist_src_ptr, &m_dists_sources[0], dist_src_size,
                         cudaMemcpyHostToDevice));

    return Multiplane(lens_ptr, data.size(), src_ptr, m_source_z.size(),
                      dist_lens_ptr, dist_src_ptr, m_dist_offsets, subl_size,
                      true);
}

Multiplane *MultiplaneBuilder::getCuMultiPlanePtr() {
    std::vector<CompositeLens> data;
    std::vector<int> subl_size;

    for (size_t i = 0; i < m_builders.size(); i++) {
        subl_size.push_back(m_builders[i].length());
        auto lens = m_builders[i].getCuLens();
        data.push_back(lens);
    }

    if (data.size() == 0 || m_source_z.size() == 0) {
        std::cerr << "No lens and/or source planes given " << data.size() << "-"
                  << m_source_z.size() << std::endl;
        throw(-1);
    }

    prepare();

    CompositeLens *lens_ptr;
    float *src_ptr, *dist_lens_ptr, *dist_src_ptr;

    size_t lens_size = sizeof(CompositeLens) * data.size();
    gpuErrchk(cudaMalloc(&lens_ptr, lens_size));
    gpuErrchk(
        cudaMemcpy(lens_ptr, &data[0], lens_size, cudaMemcpyHostToDevice));

    size_t src_size = sizeof(float) * m_source_z.size();
    gpuErrchk(cudaMalloc(&src_ptr, src_size));
    gpuErrchk(
        cudaMemcpy(src_ptr, &m_source_z[0], src_size, cudaMemcpyHostToDevice));

    size_t dist_lens_size = sizeof(float) * m_dists_lenses.size();
    gpuErrchk(cudaMalloc(&dist_lens_ptr, dist_lens_size));
    gpuErrchk(cudaMemcpy(dist_lens_ptr, &m_dists_lenses[0], dist_lens_size,
                         cudaMemcpyHostToDevice));

    size_t dist_src_size = sizeof(float) * m_dists_sources.size();
    gpuErrchk(cudaMalloc(&dist_src_ptr, dist_src_size));
    gpuErrchk(cudaMemcpy(dist_src_ptr, &m_dists_sources[0], dist_src_size,
                         cudaMemcpyHostToDevice));

    return new Multiplane(lens_ptr, data.size(), src_ptr, m_source_z.size(),
                          dist_lens_ptr, dist_src_ptr, m_dist_offsets,
                          subl_size, true);
}

int Multiplane::destroy() {
    // Destroy children
    if (m_cuda) {
        // So we copy them back to cpu first and then destroy TODO:
        // consider just keeping a vector of them in memory for
        // cleanliness
        size_t psize = m_lenses_size * sizeof(CompositeLens);
        CompositeLens *pptr = (CompositeLens *)malloc(psize);
        cpuErrchk(pptr);
        gpuErrchk(cudaMemcpy(pptr, m_lenses, psize, cudaMemcpyDeviceToHost));
        for (int i = 0; i < m_lenses_size; i++) {
            pptr[i].destroy();
        }
        free(pptr);
    } else {
        for (int i = 0; i < m_lenses_size; i++) {
            m_lenses[i].destroy();
        }
    }

    // Free memory
    if (m_cuda) {
        gpuErrchk(cudaFree(m_lenses));
        gpuErrchk(cudaFree(m_sources));
        gpuErrchk(cudaFree(m_dist_lenses));
        gpuErrchk(cudaFree(m_dist_sources));
        if (m_alphas)
            gpuErrchk(cudaFree(m_alphas));
    } else {
        free(m_lenses);
        free(m_sources);
        free(m_dist_lenses);
        free(m_dist_sources);
    }
    m_lenses = nullptr;
    m_sources = nullptr;
    m_dist_lenses = nullptr;
    m_dist_sources = nullptr;
    m_alphas = nullptr;
    return 0;
}

__global__ void mp_traceThetaGlobal(const long n,
                                    const float2 *__restrict__ thetas,
                                    float2 *__restrict__ betas,
                                    const float *__restrict__ dist_lenses,
                                    const float *__restrict__ dist_sources,
                                    const int numlenses, const int offset,
                                    const CompositeLens *__restrict__ lenses) {
    const long z = min((long)blockIdx.x * blockDim.x + threadIdx.x, n - 1);

    const int MAX_PLANES = 128;
    float2 alphas[MAX_PLANES];

    float2 last_theta;
    int l = 0;
    for (int i = 0; i <= numlenses; i++) {
        auto t = thetas[z];
        if (i > 0) {
            const int idx = i - 1;
            alphas[idx] = lenses[idx].getAlpha(last_theta);
        }
        for (int j = 0; j < i; j++) {
            const float2 a = alphas[j];
            const float dist = dist_lenses[l];
            t.x -= a.x * dist;
            t.y -= a.y * dist;
            l++;
        }
        last_theta = t;
    }

    l = offset;
    auto t = thetas[z];
    for (int i = 0; i < numlenses; i++) {
        const float2 a = alphas[i];
        const float dist = dist_sources[l];
        t.x -= a.x * dist;
        t.y -= a.y * dist;
        l++;
    }
    betas[z] = t;
}

int Multiplane::traceThetas(const float2 *thetas, float2 *betas, const size_t n,
                            const int plane) const {
    size_t offset = 0;
    for (int i = 0; i < plane; i++) {
        int s = m_dist_offsets[i];
        offset += s;
    }

    int numlenses = m_dist_offsets[plane];
    mp_traceThetaGlobal<<<(n / 128) + 1, 128>>>(n, thetas, betas, m_dist_lenses,
                                                m_dist_sources, numlenses,
                                                offset, m_lenses);
    gpuErrchk(cudaGetLastError());

    return 0;
}

__global__ void mp_updateMasses(const size_t n, const float *__restrict__ masses,
                                const int lens,
                                CompositeLens *__restrict__ lenses) {
    const size_t i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < n) {
        // printf("Mass %i: %f\n", i, masses[i]);
        lenses[lens].update(masses[i], i);
    }
}

void Multiplane::updateMassesCu(const std::vector<std::vector<float>> &masses) {
    thrust::device_vector<float> mass;
    float *ptr;
    for (size_t i = 0; i < masses.size(); i++) {
        mass = masses[i];
        ptr = thrust::raw_pointer_cast(&mass[0]);
        size_t size = masses[i].size();
        mp_updateMasses<<<(size / 64) + 1, 64>>>(size, ptr, i, m_lenses);
    }
}
