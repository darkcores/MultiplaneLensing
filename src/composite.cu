#include "composite.h"
#include "util/error.h"
#include <iostream>

int CompositeLens::destroy() {
    if (m_data_ptr) {
        if (m_cuda) {
            gpuErrchk(cudaFree(m_data_ptr));
        } else {
            free(m_data_ptr);
        }
        m_data_ptr = NULL;
    }
    return 0;
}

CompositeLens CompositeLensBuilder::getCuLens() {
#ifdef __CUDACC__
    // m_lenses.push_back(LensData());
    // dev_m_lenses = m_lenses;
    // auto lens_ptr = thrust::raw_pointer_cast(&dev_m_lenses[0]);
    size_t size = m_lenses.size();
    if (size == 0) {
        std::cerr << "No lenses added" << std::endl;
        std::terminate();
    }
    gpuErrchk(cudaMalloc(&lens_ptr, sizeof(LensData) * size));
    gpuErrchk(cudaMemcpy(lens_ptr, &m_lenses[0], sizeof(LensData) * size,
                         cudaMemcpyHostToDevice));
#else
    CompositeLens *lens_ptr = nullptr;
#endif
    return CompositeLens(m_Dd, m_Ds, m_Dds, lens_ptr, size, m_scale, true);
}

#ifdef NOTHING
__global__ void sumalphas(const int n, const float2 theta, float2 sum,
                          const LensData *__restrict__ cur_data_ptr,
                          const float m_scale) {
    __shared__ float2 blocksum;
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        if (threadIdx.x == 0) {
            blocksum.x = 0;
            blocksum.y = 0;
        }
        LensData ld = cur_data_ptr[i];
        float2 movedtheta = theta;
        movedtheta.x *= m_scale;
        movedtheta.y *= m_scale;
        movedtheta.x -= ld.position.x;
        movedtheta.y -= ld.position.y;
        movedtheta = ld.lens.getAlphaf(movedtheta);
        atomicAdd(&blocksum.x, movedtheta.x);
        atomicAdd(&blocksum.y, movedtheta.y);
        if (threadIdx.x == 0) {
            // atomicAdd(&(sum->x), blocksum.x);
            // atomicAdd(&(sum->y), blocksum.y);
        }
    }
}

float2 CompositeLens::getAlphaf(const float2 &theta) const {
    float2 alpha, movedtheta;
    LensData ld;
    alpha.x = 0;
    alpha.y = 0;
    /*
#ifdef __CUDA_ARCH__
sumalphas<<<(length / 64) + 1, 64>>>(length, theta, alpha,
                                 cur_data_ptr, m_scale);
#endif
    */
#ifdef __CUDA_ARCH__
#pragma unroll 16
#endif
    for (int i = 0; i < length; i++) {
        ld = cur_data_ptr[i];
        movedtheta = theta;
        movedtheta.x *= m_scale;
        movedtheta.y *= m_scale;
        movedtheta.x -= ld.position.x;
        movedtheta.y -= ld.position.y;
        movedtheta = ld.lens.getAlphaf(movedtheta);
        alpha.x += movedtheta.x;
        alpha.y += movedtheta.y;
    }
    // theta /= m_scale;
    alpha.x /= m_scale;
    alpha.y /= m_scale;
    return alpha;
}
#endif
