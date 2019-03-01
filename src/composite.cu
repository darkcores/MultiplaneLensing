#include "composite.h"

CompositeLensBuilder::CompositeLensBuilder() {}

void CompositeLensBuilder::addLens(Plummer &lens, Vector2D<float> position) {
    m_lenses.push_back(LensData(lens, position));
}

void CompositeLensBuilder::clear() { m_lenses.clear(); }

void CompositeLensBuilder::setDistance(const double Dd) {
    m_Dd = Dd;
    for (size_t i = 0; i < m_lenses.size(); i++) {
        m_lenses[i].lens.setDistance(Dd);
    }
}

void CompositeLensBuilder::setSource(const double Ds, const double Dds) {
    m_Ds = Ds;
    m_Dds = Dds;
    for (size_t i = 0; i < m_lenses.size(); i++) {
        m_lenses[i].lens.setSource(Ds, Dds);
    }
}

void CompositeLensBuilder::setScale(const float scale) {
    m_scale = scale;
    for (size_t i = 0; i < m_lenses.size(); i++) {
        m_lenses[i].lens.setScale(scale);
    }
}

CompositeLens CompositeLensBuilder::getLens() {
    CompositeLens lens(m_Dd, m_Ds, m_Dds, &m_lenses[0], m_lenses.size());
    return lens;
}

CompositeLens CompositeLensBuilder::getCuLens() {
#ifdef __CUDACC__
    dev_m_lenses = m_lenses;
    auto lens_ptr = thrust::raw_pointer_cast(&dev_m_lenses[0]);
#else
    CompositeLens *lens_ptr = 0;
#endif
    CompositeLens lens(m_Dd, m_Ds, m_Dds, lens_ptr, m_lenses.size(), m_scale);
    return lens;
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

__host__ __device__ Vector2D<float>
CompositeLens::getAlphaf(Vector2D<float> theta) const {
    theta *= m_scale;
    Vector2D<float> alpha(0, 0);
#ifdef __CUDA_ARCH____
    __shared__ LensData data[length];
    if (threadIdx.x == 0) {
            for (int i = 0; i < length; i++)
                    data = cur_data_ptr[i];
    }
    __syncthreads();
    for (int i = 0; i < length; i++) {
    auto movedtheta = theta - (data[i].position * m_scale);
    alpha += data[i].lens.getAlphaf(movedtheta);
    }
    #else
    for (size_t i = 0; i < length; i++) {
        /*
        #ifdef __CUDA_ARCH____
        __shared__ LensData data;
        if (threadIdx.x == 0) {
                data = cur_data_ptr[i];
        }
        __syncthreads();
auto movedtheta = theta - (data.position * m_scale);
alpha += data.lens.getAlphaf(movedtheta);
        #else
        */
        auto movedtheta = theta - (cur_data_ptr[i].position * m_scale);
        alpha += cur_data_ptr[i].lens.getAlphaf(movedtheta);
        // #endif
    }
    #endif
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
