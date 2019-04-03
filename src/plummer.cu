#include "plummer.h"

#include "util/constants.h"
#include <cmath>

#ifdef __CUDACC__
Plummer::Plummer(const double Dd, const double mass, const double angularwidth,
                 const double scale, const float2 position)
    : m_angularwidth2(angularwidth * angularwidth),
      m_4GM((4 * CONST_G * mass) / (SPEED_C * SPEED_C * Dd) * (scale * scale)),
      m_position(position) {
    m_4GM_f = m_4GM;
}
#else
Plummer::Plummer(const double Dd, const double mass, const double angularwidth,
                 const double scale, const Vector2D<float> position)
    : m_angularwidth2(angularwidth * angularwidth),
      m_4GM((4 * CONST_G * mass) / (SPEED_C * SPEED_C * Dd) * (scale * scale)),
      m_position(position) {
    m_4GM_f = m_4GM;
}
#endif
