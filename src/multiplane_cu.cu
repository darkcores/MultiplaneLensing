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
    if (m_multiplane)
        gpuErrchk(cudaFree(m_multiplane));
}

CompositeLensBuilder
MultiPlaneContext::buildLens(const double Dd,
                             const std::vector<PlummerParams> &params) {
	CompositeLensBuilder builder;
	for (auto &param : params) {
		auto position = param.position * m_angularUnit;
		Plummer plum(Dd, param.mass, param.angularwidth);
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

    return 0;
}
