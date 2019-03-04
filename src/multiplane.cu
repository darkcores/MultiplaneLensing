#include "multiplane.h"

#include <algorithm>

Multiplane::Multiplane() {}

/*
__host__ __device__ Vector2D<double> MultiPlane::getBeta(const Vector2D<double>
&theta) const {
}
*/

MultiplaneBuilder::MultiplaneBuilder() {}

void MultiplaneBuilder::addPlane(const CompositeLens &lens, double z) {
    m_data.push_back(PlaneData(lens, z));
}

void MultiplaneBuilder::addSourcePlane(const SourcePlane &plane) {
    m_src_data.push_back(plane);
}

void MultiplaneBuilder::prepare() {
    std::sort(m_data.begin(), m_data.end());
    std::sort(m_src_data.begin(), m_src_data.end());
    // Set inter lens source distances, nice that we have them sorted now
	if (m_data.size() > 1) {
		// m_data[0].setSource(
		for (size_t i = 1; i < m_data.size() - 1; i++) {
		}
	}
}
