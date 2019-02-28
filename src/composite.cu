#include "composite.h"

CompositeLens::CompositeLens() {}

void CompositeLens::addLens(Plummer &lens, Vector2D<float> position) {
    m_lenses.push_back(LensData(lens, position));
}

void CompositeLens::prepare() {
    dev_m_lenses = m_lenses;
    cur_data_ptr = thrust::raw_pointer_cast(&dev_m_lenses[0]);
    length = m_lenses.size();
}

void CompositeLens::clear() { m_lenses.clear(); }

__host__ __device__ Vector2D<double>
CompositeLens::getAlpha(Vector2D<double> theta) const {
    Vector2D<double> alpha(0, 0);
#ifdef __CUDA_ARCH__
    for (size_t i = 0; i < length; i++) {
        auto movedtheta = theta - cur_data_ptr[i].position;
        alpha += cur_data_ptr[i].lens.getAlpha(movedtheta);
    }
#else
    for (const auto &lensdata : m_lenses) {
        alpha += lensdata.lens.getAlpha(theta - lensdata.position);
    }
#endif
    return alpha;
}

__host__ __device__ Vector2D<double>
CompositeLens::getBeta(Vector2D<double> theta) const {
    Vector2D<double> beta;
    // beta = theta - (m_Dds / m_Ds) * getAlpha(theta);
	return beta;
}
