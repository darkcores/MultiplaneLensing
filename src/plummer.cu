#include "plummer.h"

#include "util/constants.h"
#include <cmath>

__host__ __device__ Plummer::Plummer(const double Dd, const double mass,
                                     const double angularwidth) {

    m_mass = mass;
    m_angularwidth = angularwidth;
    m_angularwidth2 = angularwidth * angularwidth;
    m_Dd = Dd;
    m_Ds = 0;
    m_Dds = 0;
    m_scale = 0;
    // m_angularpos = angularposition;
    m_4GM = (4 * CONST_G * m_mass) / (SPEED_C * SPEED_C * m_Dd);
    setScale();
}

__host__ __device__ Vector2D<double>
Plummer::getAlpha(Vector2D<double> theta) const {
    auto alpha = theta / (theta.lengthSq() + m_angularwidth2);
    return alpha * m_4GM;
}

__host__ __device__ Vector2D<double>
Plummer::getBeta(Vector2D<double> theta) const {
    auto beta = theta - (getAlpha(theta) * m_4GM_s);
    return beta;
}

__host__ __device__ Vector2D<float>
Plummer::getAlphaf(Vector2D<float> theta) const {
	// theta * m_scale;
    auto alpha = theta / (theta.lengthSq() + m_angularwidth2_f);
    alpha *= m_4GM_f;
	// alpha /= m_scale;
	return alpha;
}

__host__ __device__ Vector2D<float>
Plummer::getBetaf(Vector2D<float> theta) const {
    // theta *= m_scale;
    auto beta = theta / (theta.lengthSq() + m_angularwidth2_f);
    beta *= m_4GM_s_f;
    beta = theta - beta;
    // beta /= m_scale;
    return beta;
}

__host__ __device__ void Plummer::setDistance(const double Dd) {
    m_Dd = Dd;
    m_4GM = (4 * CONST_G * m_mass) / (SPEED_C * SPEED_C * m_Dd);
    m_4GM_s = (m_Dds / m_Ds) * m_4GM;
    setSource(m_Dd, m_Dds);
}

__host__ __device__ void Plummer::setSource(const double Ds, const double Dds) {
    m_Ds = Ds;
    m_Dds = Dds;
    m_4GM_s = (Dds / Ds) * m_4GM;
    setScale(m_scale);
}

__host__ __device__ void Plummer::setScale(const float scale) {
    m_scale = scale;
	m_angularwidth2_f = (m_angularwidth * m_scale) * (m_angularwidth * m_scale);
    m_4GM_f = m_4GM * (scale * scale);
    m_4GM_s_f = m_4GM_s * (scale * scale);
}

__host__ __device__ void Plummer::setMass(const double mass) {
    m_mass = mass;
    // TODO if this is slow optimize calculations
    setDistance(m_Dd);
}
