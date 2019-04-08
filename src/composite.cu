#include "composite.h"
#include "util/error.h"
#include <iostream>

__global__ void cmp_init_lenses(const int n, const Plummer *__restrict__ p,
                                float4 *__restrict__ d) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        d[i] = p[i].m_data;
    }
}

CompositeLens::CompositeLens(Plummer *lenses, const int lenses_size,
                             const bool cuda)
    : m_lenses(lenses), m_lenses_size(lenses_size), m_cuda(cuda) {
    if (m_cuda) {
        float4 *tmp;
        gpuErrchk(cudaMalloc(&tmp, sizeof(float4) * lenses_size));
        m_lens_int = tmp;
        cmp_init_lenses<<<(lenses_size / 256) + 1, 256>>>(lenses_size, m_lenses,
                                                          m_lens_int);
    }
}

int CompositeLens::destroy() {
    if (m_cuda) {
        gpuErrchk(cudaFree(m_lenses));
        gpuErrchk(cudaFree(m_lens_int));
    } else {
        free(m_lenses);
    }
    return 0;
}

CompositeLens CompositeLensBuilder::getCuLens() {
    // printf("Composite size: %lu\n", sizeof(CompositeLens));
    // printf("Plummer size: %lu\n", sizeof(Plummer));
    Plummer *lens_ptr;
    size_t size = m_lenses.size();
    // printf("Lenses in composite: %lu\n", size);
    if (size == 0) {
        std::cerr << "No lenses added" << std::endl;
        throw(-1);
    }
    size_t numbytes = sizeof(Plummer) * size;
    // printf("Bytes: %lu\n", numbytes);
    gpuErrchk(cudaMalloc(&lens_ptr, numbytes));
    gpuErrchk(
        cudaMemcpy(lens_ptr, &m_lenses[0], numbytes, cudaMemcpyHostToDevice));
    return CompositeLens(lens_ptr, (int)size, true);
}
