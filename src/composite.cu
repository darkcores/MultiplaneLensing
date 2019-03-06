#include "composite.h"
#include <iostream>

CompositeLens CompositeLensBuilder::getCuLens() {
	cuda = true;
#ifdef __CUDACC__
    // m_lenses.push_back(LensData());
    // dev_m_lenses = m_lenses;
    // auto lens_ptr = thrust::raw_pointer_cast(&dev_m_lenses[0]);
    size_t size = m_lenses.size();
    if (cudaMalloc(&lens_ptr, sizeof(LensData) * size) != cudaSuccess) {
		std::cout << "Malloc: " << cudaGetErrorString(cudaPeekAtLastError()) << std::endl;
        std::terminate();
    }
    if (cudaMemcpy(lens_ptr, &m_lenses[0], sizeof(LensData) * size,
                   cudaMemcpyHostToDevice) != cudaSuccess) {
		std::cout << "Memcpy: " << cudaGetErrorString(cudaPeekAtLastError()) << std::endl;
        std::terminate();
    }
#else
    CompositeLens *lens_ptr = 0;
#endif
    CompositeLens lens(m_Dd, m_Ds, m_Dds, lens_ptr, m_lenses.size(), m_scale);
    return lens;
}

void CompositeLensBuilder::cuFree() {
	cudaFree(lens_ptr);
}

CompositeLens::CompositeLens(const double Dd, const double Ds, const double Dds,
                             LensData *data_ptr, size_t size, float scale)
    : cur_data_ptr(data_ptr) {
    m_Dd = Dd;
    m_Ds = Ds;
    m_Dds = Dds;
    m_D = m_Dds / m_Ds;
    m_Df = m_Dds / m_Ds;
    m_scale = scale;
    // cur_data_ptr = data_ptr;
    m_data_ptr = data_ptr;
    length = size;
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