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

__host__ __device__ Vector2D<double>
CompositeLens::getAlpha(Vector2D<double> theta) const {
    Vector2D<double> alpha(0, 0);
    for (size_t i = 0; i < length; i++) {
        auto movedtheta = theta - cur_data_ptr[i].position;
        alpha += cur_data_ptr[i].lens.getAlpha(movedtheta);
    }
    return alpha;
}

__host__ __device__ Vector2D<double>
CompositeLens::getBeta(Vector2D<double> theta) const {
    Vector2D<double> beta;
    beta = theta - getAlpha(theta) * m_D;
    return beta;
}

/*
__host__ __device__ Vector2D<float>
CompositeLens::getAlphaf(Vector2D<float> theta) const {
    theta *= m_scale;
    Vector2D<float> alpha(0, 0);
    for (size_t i = 0; i < length; i++) {
        auto movedtheta = theta - (cur_data_ptr[i].position * m_scale);
        alpha += cur_data_ptr[i].lens.getAlphaf(movedtheta);
    }
    // theta /= m_scale;
    alpha /= m_scale;
    return alpha;
}

__host__ __device__ Vector2D<float>
CompositeLens::getBetaf(Vector2D<float> theta) const {
    Vector2D<float> beta;
    beta = theta - getAlphaf(theta) * m_Df;
    return beta;
}
*/