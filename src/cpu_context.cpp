#include "cpu_context.h"

#include "multiplane.h"
#include <algorithm>
#include <chrono>
#include <cstring>
#include <numeric>

CPUMultiPlaneContext::CPUMultiPlaneContext(const double angularUnit,
                                           const Cosmology cosmology)
    : m_angularUnit(angularUnit), m_cosmology(cosmology) {
    m_theta = nullptr;
    m_beta = nullptr;
    m_multiplane = nullptr;
}

CPUMultiPlaneContext::~CPUMultiPlaneContext() {
    if (m_theta)
        free(m_theta);
    if (m_beta)
        free(m_beta);
    if (m_multiplane) {
        m_multiplane->destroy();
        free(m_multiplane);
    }
}

CompositeLensBuilder
CPUMultiPlaneContext::buildLens(const float redshift,
                                const std::vector<PlummerParams> &params) {
    CompositeLensBuilder builder(redshift);
    double Dd = m_cosmology.angularDiameterDistance(redshift);
    for (auto &param : params) {
        Plummer plum(Dd, param.mass, param.angularwidth, 1 / m_angularUnit,
                     param.position);
        builder.addLens(plum);
    }
    return builder;
}

int CPUMultiPlaneContext::init(
    const std::vector<float> &lensRedshifts,
    const std::vector<std::vector<PlummerParams>> &params,
    const std::vector<float> &sourceRedshifts) {
    try {
        MultiplaneBuilder planebuilder(m_cosmology);

        // Setup lenses
        for (size_t i = 0; i < lensRedshifts.size(); i++) {
            auto lens = buildLens(lensRedshifts[i], params[i]);
            planebuilder.addPlane(lens);
        }

        // Setup sources
        planebuilder.setRedshifts(sourceRedshifts);

        // Build multiplane
        m_multiplane = planebuilder.getMultiPlanePtr();

        return 0;
    } catch (int e) {
        return e;
    }
}

int CPUMultiPlaneContext::setThetas(
    const std::vector<std::vector<Vector2D<float>>> &thetas) {
    try {
        m_theta_len = 0;
        // Resize list of betas
        std::cout << "Thetas size " << thetas.size() << std::endl;
        m_betas.resize(thetas.size());

        for (size_t i = 0; i < thetas.size(); i++) {
            size_t s = thetas[i].size();
            // printf("Theta (%lu) len: %lu\n", i, s);
            m_theta_len += s;
            m_theta_count.push_back(m_theta_len);
            m_betas[i].reserve(s);
        }

        // printf("Theta len: %lu\n", m_theta_len);

        m_theta = (Vector2D<float> *)std::malloc(sizeof(Vector2D<float>) *
                                                 m_theta_len);
        size_t offset = 0;
        for (size_t i = 0; i < thetas.size(); i++) {
            size_t size = sizeof(Vector2D<float>) * thetas[i].size();
            std::memcpy(&m_theta[offset], &thetas[i][0], size);
            offset = m_theta_count[i];
        }
        size_t beta_size = sizeof(Vector2D<float>) * m_theta_len;
        m_beta = (Vector2D<float> *)malloc(beta_size);

        return 0;
    } catch (int e) {
        return e;
    }
}

int CPUMultiPlaneContext::calculatePositionsBenchmark(
    const std::vector<std::vector<float>> &masses, float &millis, int nruns) {
    try {
        std::vector<float> m(nruns);

        for (int x = 0; x < nruns; x++) {

            // Setup new masses
            m_multiplane->updateMasses(masses);

            // Calculate new betas
            auto start = std::chrono::system_clock::now();
            size_t offset = 0;
            for (size_t i = 0; i < m_betas.size(); i++) {
                size_t tcount = m_theta_count[i] - offset;
                m_multiplane->traceThetas((Vector2D<float> *)&m_theta[offset],
                                          (Vector2D<float> *)&m_beta[offset],
                                          tcount, i);

                // Copy results back to host
                m_betas[i].resize(tcount);
                std::memcpy(&m_betas[i][0], &m_beta[offset],
                            sizeof(Vector2D<float>) * tcount);
                offset = m_theta_count[i];
            }
            auto end = std::chrono::system_clock::now();

            std::chrono::duration<double> elapsed_seconds = end - start;
            elapsed_seconds *= 1000;
            m[x] = elapsed_seconds.count();
            // m[x] =
            // std::chrono::duration_cast<std::chrono::milliseconds>(elapsed_seconds).count();
        }

        float avg = std::accumulate(m.begin(), m.end(), 0) / (float)m.size();
        millis = avg;

        float minv = *std::min_element(m.begin(), m.end());
        float maxv = *std::max_element(m.begin(), m.end());
        float med = m[m.size() / 2];

        float sq_sum = std::inner_product(m.begin(), m.end(), m.begin(), 0.0);
        float stdev = std::sqrt(sq_sum / m.size() - avg * avg);

        printf(
            "Avg time: %.2f ms (min: %.2f; med: %.2f; max: %.2f, std: %.2f)\n",
            avg, minv, med, maxv, stdev);

        return 0;
    } catch (int e) {
        return e;
    }
}

const std::vector<Vector2D<float>> &
CPUMultiPlaneContext::getSourcePositions(int idx) const {
    return m_betas[idx];
}
