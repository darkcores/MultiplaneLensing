#include "plummer.h"

#include "util/constants.h"
#include <cmath>

__host__ __device__ Plummer::Plummer(const double Dd, const double mass,
                                     const double angularwidth) {

    m_mass = mass;
    m_angularwidth = angularwidth;
    m_Dd = Dd;
	m_Df = 1;
    m_scale = 1;
    // m_angularpos = angularposition;
    setScale(1);
}

__host__ __device__ void Plummer::setDistance(const double Dd) {
    m_Dd = Dd;
    setScale(m_scale);
}

__host__ __device__ void Plummer::setSource(const double Ds, const double Dds) {
	m_Df = Dds / Ds;
    setScale(m_scale);
}

__host__ __device__ void Plummer::setScale(const float scale) {
    m_scale = scale;
    m_angularwidth2_f = (m_angularwidth * m_scale) * (m_angularwidth * m_scale);
    m_4GM_f = (4 * CONST_G * m_mass) / (SPEED_C * SPEED_C * m_Dd) * (scale * scale);
}

__host__ __device__ void Plummer::setMass(const double mass) {
    m_mass = mass;
    // printf("Setting mass on lens: %f\n", mass);
    // TODO if this is slow optimize calculations
    setDistance(m_Dd);
}
