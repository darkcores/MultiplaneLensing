#include "multiplane_cu.h"

#include "multiplane.h"
#include "util/error.h"

MultiPlaneContext::MultiPlaneContext(const double angularUnit,
                                     const Cosmology cosmology)
    : m_angularUnit(angularUnit), m_cosmology(cosmology) {
    m_theta_x = nullptr;
    m_theta_y = nullptr;
    m_beta_x = nullptr;
    m_beta_y = nullptr;
    m_multiplane = nullptr;
}

MultiPlaneContext::~MultiPlaneContext() {
    if (m_theta_x)
        gpuErrchk(cudaFree(m_theta_x));
    if (m_theta_y)
        gpuErrchk(cudaFree(m_theta_y));
    if (m_beta_x)
        gpuErrchk(cudaFree(m_beta_x));
    if (m_beta_y)
        gpuErrchk(cudaFree(m_beta_y));
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
    MultiplaneBuilder planebuilder(m_cosmology);

    // Setup lenses
    for (size_t i = 0; i < lensRedshifts.size(); i++) {
        double Dd = m_cosmology.angularDiameterDistance(lensRedshifts[i]);
        auto lens = buildLens(Dd, params[i]);
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

    return 0;
}

int MultiPlaneContext::setThetas(const std::vector<Vector2D<float>> &thetas) {
    // Temp cpu arrays
    std::vector<float> theta_x, theta_y;
    // Allocations for cuda
    size_t arr_size = sizeof(float) * thetas.size();
    gpuErrchk(cudaMalloc(&m_theta_x, arr_size));
    gpuErrchk(cudaMalloc(&m_theta_y, arr_size));
    size_t beta_size = arr_size * m_multiplane->srcLength();
    gpuErrchk(cudaMalloc(&m_beta_x, beta_size));
    gpuErrchk(cudaMalloc(&m_beta_y, beta_size));

    for (auto theta : thetas) {
		theta *= m_angularUnit;
        theta_x.push_back(theta.x());
        theta_y.push_back(theta.y());
    }

    gpuErrchk(
        cudaMemcpy(m_theta_x, &theta_x[0], arr_size, cudaMemcpyHostToDevice));
    gpuErrchk(
        cudaMemcpy(m_theta_y, &theta_y[0], arr_size, cudaMemcpyHostToDevice));

    return 0;
}
