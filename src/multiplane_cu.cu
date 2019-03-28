#include "multiplane_cu.h"

#include "multiplane.h"
#include "util/error.h"

MultiPlaneContext::MultiPlaneContext(const double angularUnit,
                                     const Cosmology cosmology)
    : m_angularUnit(angularUnit), m_cosmology(cosmology) {
    // m_theta_x = nullptr;
    // m_theta_y = nullptr;
    m_theta = nullptr;
    // m_beta_x = nullptr;
    // m_beta_y = nullptr;
    m_beta = nullptr;
    m_multiplane = nullptr;

	// gpuErrchk(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1));
}

MultiPlaneContext::~MultiPlaneContext() {
    /*
if (m_theta_x)
    gpuErrchk(cudaFree(m_theta_x));
if (m_theta_y)
    gpuErrchk(cudaFree(m_theta_y));
    */
    if (m_theta)
        gpuErrchk(cudaFree(m_theta));
    /*
if (m_beta_x)
gpuErrchk(cudaFree(m_beta_x));
if (m_beta_y)
gpuErrchk(cudaFree(m_beta_y));
    */
    if (m_beta)
        gpuErrchk(cudaFree(m_beta));
    if (m_multiplane) {
        m_multiplane->destroy();
        free(m_multiplane);
    }
}

CompositeLensBuilder
MultiPlaneContext::buildLens(const double Dd,
                             const std::vector<PlummerParams> &params) {
    CompositeLensBuilder builder;
    for (auto &param : params) {
        auto position = param.position * m_angularUnit;
        Plummer plum(Dd, param.mass, param.angularwidth * m_angularUnit);
        builder.addLens(plum, position);
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
            double Dd = m_cosmology.angularDiameterDistance(lensRedshifts[i]);
            auto lens = buildLens(Dd, params[i]);
            lens.setScale(3600);
            lens.setRedshift(lensRedshifts[i]);
            planebuilder.addPlane(lens);
        }

        // Setup sources
        for (auto z : sourceRedshifts) {
            // TODO, maybe just add a function that takes this vector
            SourcePlaneBuilder sp(z);
            auto plane = sp.getCuPlane();
            planebuilder.addSourcePlane(plane);
        }

        // Build multiplane
        m_multiplane = (Multiplane *)malloc(sizeof(Multiplane));
        cpuErrchk(m_multiplane);
        auto plane = planebuilder.getCuMultiPlane();
        memcpy(m_multiplane, &plane, sizeof(Multiplane));

        // Resize list of betas
        m_betas.resize(sourceRedshifts.size());

        return 0;
    } catch (int e) {
        return e;
    }
}

int MultiPlaneContext::setThetas(const std::vector<Vector2D<float>> &thetas) {
    try {
        /*
// Temp cpu arrays
std::vector<float> theta_x, theta_y;
// Allocations for cuda
size_t arr_size = sizeof(float) * thetas.size();
gpuErrchk(cudaMalloc(&m_theta_x, arr_size));
gpuErrchk(cudaMalloc(&m_theta_y, arr_size));
size_t beta_size = arr_size * m_multiplane->srcLength();
gpuErrchk(cudaMalloc(&m_beta_x, beta_size));
gpuErrchk(cudaMalloc(&m_beta_y, beta_size));

m_theta_len = thetas.size();
for (auto theta : thetas) {
    theta *= m_angularUnit;
    theta_x.push_back(theta.x());
    theta_y.push_back(theta.y());
}

gpuErrchk(cudaMemcpy(m_theta_x, &theta_x[0], arr_size,
                     cudaMemcpyHostToDevice));
gpuErrchk(cudaMemcpy(m_theta_y, &theta_y[0], arr_size,
                     cudaMemcpyHostToDevice));

// Reserve space for betas
for (auto &x : m_betas) {
    x.reserve(thetas.size());
}
        */
        m_theta_len = thetas.size();
		std::vector<float2> tmp;
		tmp.reserve(m_theta_len);
		for (auto &x : thetas) {
			float2 f;
			f.x = x.x() * m_angularUnit;
			f.y = x.y() * m_angularUnit;
			tmp.push_back(f);
		}
        size_t size = sizeof(Vector2D<float>) * thetas.size();
        gpuErrchk(cudaMalloc(&m_theta, size));
        size_t beta_size = size * m_multiplane->srcLength();
        gpuErrchk(cudaMalloc(&m_beta, beta_size));
        gpuErrchk(
            cudaMemcpy(m_theta, &tmp[0], size, cudaMemcpyHostToDevice));

        // Reserve space for betas
        for (auto &x : m_betas) {
            x.reserve(thetas.size());
        }

        return 0;
    } catch (int e) {
        return e;
    }
}

/**
 * Where n is the number of masses / lenses and plane the lensplane.
 */
__global__ void _updateMassesKernel(const int n, const int plane, Multiplane mp,
                                    const double *masses) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n)
        mp.updateLensMasses(plane, i, masses);
}

__global__ void _traceThetaKernel(const int n, const Multiplane mp,
                                  const float2 *__restrict__ points,
                                  float2 *__restrict__ betas) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        // Vector2D<float> vec(xpoints[i], ypoints[i]);
        mp.traceTheta(points[i], &betas[i], n);
    }
}

int MultiPlaneContext::calculatePositions(
    const std::vector<std::vector<double>> &masses) {
    try {
        // Setup new masses
        size_t size = masses.size();
        double *mass;
        size_t lsize = 0;
        for (auto &m : masses) {
            if (m.size() > lsize)
                lsize = m.size();
        }
        // This array could be permanent if malloc / free is too long
        gpuErrchk(cudaMalloc(&mass, lsize * sizeof(double)));
        for (size_t i = 0; i < size; i++) {
            size_t msize = masses[i].size();
            gpuErrchk(cudaMemcpy(mass, &masses[i][0], sizeof(double) * msize,
                                 cudaMemcpyHostToDevice));

            _updateMassesKernel<<<(size / 32) + 1, 32>>>(msize, i,
                                                         *m_multiplane, mass);
        }
        gpuErrchk(cudaFree(mass));

        // Calculate new betas
        _traceThetaKernel<<<(m_theta_len / 256) + 1, 256>>>(
            m_theta_len, *m_multiplane, (float2 *)m_theta, (float2 *)m_beta);

        // Copy results back to host
        for (size_t i = 0; i < m_betas.size(); i++) {
            m_betas[i].resize(m_theta_len);
            gpuErrchk(cudaMemcpy(&m_betas[i][0], &m_beta[i * m_theta_len],
                                 sizeof(float2) * m_theta_len,
                                 cudaMemcpyDeviceToHost));
			for (auto &x : m_betas[i]) {
				x /= m_angularUnit;
			}
        }

        return 0;
    } catch (int e) {
        return e;
    }
}

const std::vector<Vector2D<float>> &
MultiPlaneContext::getSourcePositions(int idx) const {
    return m_betas[idx];
}
