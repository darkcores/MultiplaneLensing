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

/*
void Multiplane::alphaSetup() {
    auto max =
        std::max_element(m_dist_offsets.begin(), m_dist_offsets.end());
        int n = m_sources_size;
    float *atmp = nullptr;
    gpuErrchk(cudaMalloc(&atmp, n * sizeof(float2) * *max));
    m_alphas = atmp;
}
*/

void Multiplane::prepare(int numThetas) {
    if (m_alphas != nullptr)
        gpuErrchk(cudaFree(m_alphas));

    float *atmp = nullptr;
    int numlenses =
        *std::max_element(m_dist_offsets.begin(), m_dist_offsets.end());
    // printf("Alloc %d bytes\n", numlenses * sizeof(float2) * numThetas);
    gpuErrchk(cudaMalloc(&atmp, numThetas * sizeof(float2) * numlenses));
    m_alphas = atmp;
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

__global__ void mp_traceThetaGlobal(const int n,
                                    const float2 *__restrict__ thetas,
                                    float2 *__restrict__ betas, const int plane,
                                    const float *__restrict__ dist_lenses,
                                    const float *__restrict__ dist_sources,
                                    const int numlenses, const int offset,
                                    const CompositeLens *__restrict__ lenses,
                                    float2 *__restrict__ alphas) {
    unsigned int z = blockIdx.x * blockDim.x + threadIdx.x;

    if (z >= n) {
        z = n - 1;
    }
    float2 last_theta;
    // const int alphaoffset = z * numlenses;
    unsigned int l = 0;
    unsigned int idj;
    for (unsigned int i = 0; i <= numlenses; i++) {
        auto t = thetas[z];
        if (i > 0) {
            idj = z + (n * (i - 1));
            alphas[idj] = lenses[i - 1].getAlpha(last_theta);
        }
        idj = z;
        for (unsigned int j = 0; j < i; j++) {
            const float2 a = alphas[idj];
            const float dist = dist_lenses[l];
            t.x -= a.x * dist;
            t.y -= a.y * dist;
            l++;
            idj += n;
        }
        last_theta = t;
    }

    l = offset;
    auto t = thetas[z];
    idj = z;
    for (unsigned int i = 0; i < (numlenses); i++) {
        const float2 a = alphas[idj];
        const float dist = dist_sources[l];
        t.x -= a.x * dist;
        t.y -= a.y * dist;
        l++;
        idj += n;
    }
    betas[z] = t;
}

int Multiplane::traceThetas(const float2 *thetas, float2 *betas, const int n,
                            const int plane) const {
    int offset = 0;
    for (int i = 0; i < plane; i++) {
        int s = m_dist_offsets[i];
        offset += s;
    }

    int numlenses = m_dist_offsets[plane];
    if (m_alphas == nullptr) {
        thrust::device_vector<float2> alphas(numlenses * n);
        float2 *ptr = thrust::raw_pointer_cast(&alphas[0]);
        mp_traceThetaGlobal<<<(n / 128) + 1, 128>>>(
            n, thetas, betas, plane, m_dist_lenses, m_dist_sources, numlenses,
            offset, m_lenses, ptr);
    } else {
        // gpuErrchk(cudaMemset(m_alphas, 0, n * numlenses * sizeof(float2)));
        mp_traceThetaGlobal<<<(n / 128) + 1, 128>>>(
            n, thetas, betas, plane, m_dist_lenses, m_dist_sources, numlenses,
            offset, m_lenses, (float2 *)m_alphas);
    }
    gpuErrchk(cudaGetLastError());

    return 0;
}

__global__ void mp_updateMasses(const int n, const float *__restrict__ masses,
                                const int lens,
                                CompositeLens *__restrict__ lenses) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;

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