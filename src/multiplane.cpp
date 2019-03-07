#include "multiplane.h"

#include <algorithm>
#include <iostream>

void MultiplaneBuilder::addPlane(CompositeLensBuilder *lensbuilder) {
    m_builders.push_back(lensbuilder);
}

void MultiplaneBuilder::addSourcePlane(SourcePlane &plane) {
    m_src_data.push_back(plane);
}

void MultiplaneBuilder::prepare() {
    std::sort(m_data.begin(), m_data.end());
    std::sort(m_src_data.begin(), m_src_data.end());

    // Set lens distances and stuff
    size_t s = 0; // source plane index
    float z_src = m_src_data[s].redshift();
    std::cout << "Lenses in builder " << m_builders.size() << std::endl;
    for (size_t i = 0; i < m_builders.size() - 1; i++) {
        float zd = m_builders[i]->redshift();
        float zs = m_builders[i + 1]->redshift();
        float Ds = m_cosm.angularDiameterDistance(zs);
        float Dds = m_cosm.angularDiameterDistance(zs, zd);
        std::cout << zd << " - " << zs << " | " << Ds << " - " << Dds
                  << std::endl;
        m_builders[i]->setSource(Ds, Dds);
        // For each src plane before the next lens
        while (z_src < zs) {
            // Ds and Dds for this lens for the following sourceplane
            // std::cout << "Settings redshift source plane" << std::endl;
            if (z_src > zd) {
                Ds = m_cosm.angularDiameterDistance(z_src);
                Dds = m_cosm.angularDiameterDistance(z_src, zd);
                m_src_data[s].setSource(Ds, Dds);
            }
            if (s < (m_src_data.size() - 1)) {
                s++;
                z_src = m_src_data[s].redshift();
            } else {
                // z_src = std::numeric_limits<float>::infinity();
                // Later lenses aren't useful anyways
				return;
            }
        }
    }
    // Handle leftover source plane for last lens
    while (s < m_src_data.size()) {
        std::cout << "Settings redshift source plane" << std::endl;
        float zs = m_builders[m_builders.size() - 1]->redshift();
        z_src = m_src_data[s].redshift();
        float Ds = m_cosm.angularDiameterDistance(z_src);
        float Dds = m_cosm.angularDiameterDistance(z_src, zs);
        m_src_data[s].setSource(Ds, Dds);
        s++;
    }
}

Multiplane MultiplaneBuilder::getMultiPlane() {
    prepare();

    // Get final lenses from builders
    for (size_t i = 0; i < m_builders.size(); i++) {
        auto lens = m_builders[i]->getLens();
        PlaneData plane(lens, m_builders[i]->redshift());
        m_data.push_back(plane);
    }

    return Multiplane(m_data.size(), m_src_data.size(), &m_data[0],
                      &m_src_data[0]);
}
