#include "composite.h"
#include <iostream>

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
    CompositeLens lens(m_Dd, m_Ds, m_Dds, &m_lenses[0], m_lenses.size(), m_scale);
    return lens;
}
