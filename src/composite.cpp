#include "composite.h"
#include <iostream>
#include <cstdlib>
#include <cstring>

void CompositeLensBuilder::addLens(Plummer &lens, Vector2D<float> position) {
    lens.setScale(m_scale);
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
    // std::cout << "Lenses in composite: " << m_lenses.size() << std::endl;
    for (size_t i = 0; i < m_lenses.size(); i++) {
        m_lenses[i].lens.setSource(Ds, Dds);
    }
}

void CompositeLensBuilder::setScale(const float scale) {
    for (size_t i = 0; i < m_lenses.size(); i++) {
        m_lenses[i].lens.setScale(scale);
        // m_lenses[i].position /= m_scale;
        // m_lenses[i].position *= scale;
    }
    m_scale = scale;
}

CompositeLens CompositeLensBuilder::getLens() {
    // m_lenses.push_back(LensData());
	lens_ptr = (LensData *)malloc(sizeof(LensData) * m_lenses.size());
	if (lens_ptr == NULL) {
		std::cerr << "Unable to allocate memory" << std::endl;
		std::terminate();
	}
	std::memcpy(lens_ptr, &m_lenses[0], sizeof(LensData) * m_lenses.size());
    CompositeLens lens(m_Dd, m_Ds, m_Dds, lens_ptr, m_lenses.size(),
                       m_scale);
    return lens;
}

CompositeLens::CompositeLens(const double Dd, const double Ds, const double Dds,
                             LensData *data_ptr, size_t size, float scale,
                             bool cuda)
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
	m_cuda = cuda;
}
