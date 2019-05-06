#include "context.h"

#include "multiplane.h"
#include <algorithm>
#include <numeric>

MultiPlaneContext::MultiPlaneContext(const double angularUnit,
                                     const Cosmology cosmology)
    : m_angularUnit(angularUnit), m_cosmology(cosmology) {
    m_theta = nullptr;
    m_beta = nullptr;
    m_multiplane = nullptr;
}

MultiPlaneContext::~MultiPlaneContext() {
    if (m_theta)
        gpuErrchk(cudaFree(m_theta));
    if (m_beta)
        gpuErrchk(cudaFree(m_beta));
    if (m_multiplane) {
        m_multiplane->destroy();
        free(m_multiplane);
    }
}

CompositeLensBuilder
MultiPlaneContext::buildLens(const float redshift,
                             const std::vector<PlummerParams> &params) {
    CompositeLensBuilder builder(redshift);
    double Dd = m_cosmology.angularDiameterDistance(redshift);
    for (auto &param : params) {
        float2 position = param.position.f2();
        Plummer plum(Dd, param.mass, param.angularwidth, 1 / m_angularUnit,
                     position);
        builder.addLens(plum);
    }
    return builder;
}

int MultiPlaneContext::init(
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
        m_multiplane = planebuilder.getCuMultiPlanePtr();

        return 0;
    } catch (int e) {
        return e;
    }
}

int MultiPlaneContext::setThetas(
    const std::vector<std::vector<Vector2D<float>>> &thetas) {
    try {
        m_theta_len = 0;
        // Resize list of betas
		// std::cout << "Thetas size " << thetas.size() << std::endl;
        m_betas.resize(thetas.size());

        for (size_t i = 0; i < thetas.size(); i++) {
            size_t s = thetas[i].size();
            // printf("Theta (%lu) len: %lu\n", i, s);
            m_theta_len += s;
            m_theta_count.push_back(m_theta_len);
            m_betas[i].reserve(s);
        }

        // printf("Theta len: %lu\n", m_theta_len);

        gpuErrchk(cudaMalloc(&m_theta, sizeof(float2) * m_theta_len));
        size_t offset = 0;
        for (size_t i = 0; i < thetas.size(); i++) {
            size_t size = sizeof(float2) * thetas[i].size();
            gpuErrchk(cudaMemcpyAsync(&m_theta[offset], &thetas[i][0], size,
                                      cudaMemcpyHostToDevice));
            offset = m_theta_count[i];
        }
        size_t beta_size = sizeof(float2) * m_theta_len;
        gpuErrchk(cudaMalloc(&m_beta, beta_size));

        return 0;
    } catch (int e) {
        return e;
    }
}

int MultiPlaneContext::calculatePositions(
    const std::vector<std::vector<float>> &masses) {
    try {
        // Setup new masses
        m_multiplane->updateMassesCu(masses);

        // cudaStream_t stream1;
        // gpuErrchk(cudaStreamCreate(&stream1));

        // Calculate new betas
        size_t offset = 0;
        for (size_t i = 0; i < m_betas.size(); i++) {
            size_t tcount = m_theta_count[i] - offset;

            // printf("Running kernel: %d, offset %d\n", i, offset);
            m_multiplane->traceThetas((float2 *)&m_theta[offset],
                                      (float2 *)&m_beta[offset], tcount, i);

            // Copy results back to host
            m_betas[i].resize(tcount);
            // printf("Betas (%d) size: %lu\n", i, m_betas[i].size());
            // printf("Betas ptr: %p\n", &m_betas[i][0]);
            // printf("Tcount: %lu\n", tcount);
            gpuErrchk(cudaMemcpyAsync(&m_betas[i][0], &m_beta[offset],
                                      sizeof(float2) * tcount,
                                      cudaMemcpyDeviceToHost /*, stream1*/));
            offset = m_theta_count[i];
        }

        cudaDeviceSynchronize();
        // gpuErrchk(cudaStreamDestroy(stream1));

        return 0;
    } catch (int e) {
        return e;
    }
}

int MultiPlaneContext::calculatePositionsBenchmark(
    const std::vector<std::vector<float>> &masses, float &millis, int nruns) {
    try {
        std::vector<float> m(nruns);

		// First few rounds to warm up the card
        for (int x = -2; x < nruns; x++) {
            cudaEvent_t start, stop;
            cudaEventCreate(&start);
            cudaEventCreate(&stop);
            // Setup new masses
            m_multiplane->updateMassesCu(masses);

            // Calculate new betas
            cudaEventRecord(start);
            size_t offset = 0;
			// std::cout << "Source planes " << m_betas.size() << std::endl;
            for (size_t i = 0; i < m_betas.size(); i++) {
                size_t tcount = m_theta_count[i] - offset;
                m_multiplane->traceThetas((float2 *)&m_theta[offset],
                                          (float2 *)&m_beta[offset], tcount, i);

                // Copy results back to host
                m_betas[i].resize(tcount);
				/*
                std::cout << "Run " << x << "; Sizes: " << tcount << "; "
                          << m_betas[i].size() << "; " << m_theta_count[0]
                          << std::endl;
				*/
                gpuErrchk(cudaMemcpyAsync(
                    &m_betas[i][0], &m_beta[offset], sizeof(float2) * tcount,
                    cudaMemcpyDeviceToHost /*, stream1*/));
                offset = m_theta_count[i];
            }
            cudaEventRecord(stop);

            cudaDeviceSynchronize();
            cudaEventSynchronize(stop);

			if (x >= 0) {
				cudaEventElapsedTime(&m[x], start, stop);
			}
        }

		std::sort(m.begin(), m.end());

		printf("Times measured: [");
		for (auto x : m) {
			printf("%.2f, ", x);
		}
		printf("\n");
		
        float avg = std::accumulate(m.begin(), m.end(), 0.0);
		avg /= m.size();
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
MultiPlaneContext::getSourcePositions(int idx) const {
    return m_betas[idx];
}
