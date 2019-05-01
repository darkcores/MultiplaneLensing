#include "multiplane.h"

#include "util/error.h"
#include <algorithm>
#include <cstring>
#include <iostream>

void MultiplaneBuilder::addPlane(CompositeLensBuilder &lensbuilder) {
    m_builders.push_back(lensbuilder);
}

void MultiplaneBuilder::prepare() {
    std::sort(m_builders.begin(), m_builders.end());
    std::sort(m_source_z.begin(), m_source_z.end());

    // Calculate distances for all the lenses
    double sz, lz, Di, Dji, Dd;
    for (size_t i = 0; i <= m_builders.size(); i++) {
        lz = m_builders[i].redshift();
        Di = m_cosm.angularDiameterDistance(lz);
        for (size_t j = 0; j < i; j++) {
            Dji = m_cosm.angularDiameterDistance(m_builders[j].redshift(), lz);
            Dd = Dji / Di;
            m_dists_lenses.push_back(Dd);
        }
    }

    // And finally for the source planes
    for (size_t i = 0; i < m_source_z.size(); i++) {
        sz = m_source_z[i];
        Di = m_cosm.angularDiameterDistance(sz);
        size_t li = 0;
        while (li < m_builders.size() && m_builders[li].redshift() < sz) {
            Dji = m_cosm.angularDiameterDistance(m_builders[li].redshift(), sz);
            Dd = Dji / Di;
            m_dists_sources.push_back(Dd);
            li++;
        }
        m_dist_offsets.push_back(li);
        // printf("Li val %lu, %lu\n", i, li);
    }
}

Multiplane MultiplaneBuilder::getMultiPlane() {
    std::vector<CompositeLens> data;
	std::vector<int> subl_size;

    // Get final lenses from builders
    for (size_t i = 0; i < m_builders.size(); i++) {
		subl_size.push_back(m_builders[i].length());
        auto lens = m_builders[i].getLens();
        data.push_back(lens);
    }

    if (data.size() == 0 || m_source_z.size() == 0) {
        std::cerr << "No lens and/or source planes given " << data.size() << "-"
                  << m_source_z.size() << std::endl;
        throw(-1);
    }

    prepare();

    CompositeLens *lens_ptr;
    float *src_ptr, *dist_lens_ptr, *dist_src_ptr;

    size_t lens_size = sizeof(CompositeLens) * data.size();
    lens_ptr = (CompositeLens *)malloc(lens_size);
    cpuErrchk(lens_ptr);
    std::memcpy((void *)lens_ptr, &data[0], lens_size);

    size_t src_size = sizeof(float) * m_source_z.size();
    src_ptr = (float *)malloc(src_size);
    cpuErrchk(src_ptr);
    std::memcpy((void *)src_ptr, &m_source_z[0], src_size);

    size_t dist_lens_size = sizeof(float) * m_dists_lenses.size();
    dist_lens_ptr = (float *)malloc(dist_lens_size);
    cpuErrchk(dist_lens_ptr);
    std::memcpy((void *)dist_lens_ptr, &m_dists_lenses[0], dist_lens_size);

    size_t dist_src_size = sizeof(float) * m_dists_sources.size();
    dist_src_ptr = (float *)malloc(dist_src_size);
    cpuErrchk(dist_src_ptr);
    std::memcpy((void *)dist_src_ptr, &m_dists_sources[0], dist_src_size);

    return Multiplane(lens_ptr, data.size(), src_ptr, m_source_z.size(),
                      dist_lens_ptr, dist_src_ptr, m_dist_offsets, subl_size);
}

int Multiplane::traceThetas(const Vector2D<float> *thetas,
                            Vector2D<float> *betas, const size_t n,
                            const int plane) const {
    int offset = 0;
    for (int i = 0; i < plane; i++) {
        int s = m_dist_offsets[i];
        offset += s;
        // offset += (s * (s + 1) / 2);
    }

    int numlenses = m_dist_offsets[plane];
    std::vector<Vector2D<float>> alphas;
    alphas.resize(numlenses);
    // tmp_thetas.resize(numlenses);
    Vector2D<float> last_theta;

    // For each theta
    for (size_t z = 0; z < n; z++) {
        int l = 0;

        // printf("Lenses: %d\n", numlenses);

        // lenses
        for (int i = 0; i <= numlenses; i++) {
            auto t = thetas[z];
            for (int j = 0; j < i; j++) {
                if (j == (i - 1)) {
                    // Alpha not yet calculated
                    alphas[j] = m_lenses[j].getAlpha(last_theta);
                    // printf("Alpha %d\n", j);
                }
                t -= alphas[j] * m_dist_lenses[l];
                l++;
            }
            last_theta = t;
        }

        // Source plane
        l = offset;
        auto t = thetas[z];
        for (int i = 0; i < numlenses; i++) {
            t -= alphas[i] * m_dist_sources[l];
            l++;
        }
        betas[z] = t;
    }
    return 0;
}

void Multiplane::updateMasses(const std::vector<std::vector<float>> &masses) {
	for (size_t i =0; i < masses.size(); i++) {
		m_lenses[i].update(masses[i].data());
	}
}
