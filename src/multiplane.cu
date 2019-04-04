#include "multiplane.h"

#include "util/error.h"
#include <thrust/device_vector.h>
#include <algorithm>
#include <iostream>

Multiplane MultiplaneBuilder::getCuMultiPlane() {
    std::vector<CompositeLens> data;

    for (size_t i = 0; i < m_builders.size(); i++) {
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
                      dist_lens_ptr, dist_src_ptr, m_dist_offsets, true);
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
    return 0;
}

int Multiplane::traceThetas(const float2 *thetas, float2 *betas,
                            const int n, const int plane) {
    int offset = 0;
    for (int i = 0; i < plane; i++) {
        int s = m_dist_offsets[i];
		offset += s;
    }

    int numlenses = m_dist_offsets[plane];
	thrust::device_vector<float2> alphas(numlenses * n);

    return 0;
}
