#pragma once

#include "util/constants.h"
#include "util/vector2d.h"
#include <cstdio>

/**
 * Plummer lens class.
 */
class Plummer {
  private:
    const float m_angularwidth2;
    const float m_4GM;
    float m_4GM_f;
#ifdef __CUDACC__
    const float2 m_position;
#else
    const Vector2D<float> m_position;
	char __padding[4];
#endif

  public:
    /**
     * Create a new plummer lens.
     *
     * @param Dd angular distance for this lens.
     * @param mass lens mass.
     * @param angularwidth angular width for this lens.
     */
#ifdef __CUDACC__
    Plummer(const double Dd, const double mass, const double angularwidth,
            const double scale, const float2 position)
        : m_angularwidth2(angularwidth * angularwidth),
          m_4GM((4 * CONST_G * mass) / (SPEED_C * SPEED_C * Dd) *
                (scale * scale)),
          m_position(position) {
        m_4GM_f = m_4GM;
		// printf("\x1B[34m4GM: %f (m %lf)\x1B[0m \n", m_4GM_f, mass);
    }
#else
    Plummer(const double Dd, const double mass, const double angularwidth,
            const double scale, const Vector2D<float> position)
        : m_angularwidth2(angularwidth * angularwidth),
          m_4GM((4 * CONST_G * mass) / (SPEED_C * SPEED_C * Dd) *
                (scale * scale)),
          m_position(position) {
        m_4GM_f = m_4GM;
		// printf("\x1B[35m4GM: %f (m %lf)\x1B[0m \n", m_4GM_f, mass);
    }
#endif

    __host__ __device__ void update(const float scalar) {
        m_4GM_f = m_4GM * scalar;
		// printf("\x1B[31m4GM: %f (* %f)\x1B[0m \n", m_4GM_f, scalar);
    }

#ifdef __CUDACC__
    /**
     * Get alpha vector (single precision, with scaling). For scaling
     * see setScale().
     *
     * @param theta Theta vector.
     * @returns Alpha vector.
     */
    __host__ __device__ float2 getAlpha(const float2 &theta) const {
        // printf("Scale factor plummer %f\n", m_scale);
        auto alpha = theta;
        alpha.x -= m_position.x;
        alpha.y -= m_position.y;
        float len = (alpha.x * alpha.x) + (alpha.y * alpha.y);
        len += m_angularwidth2;
        len = (1 / len) * m_4GM_f;
        alpha.x *= len;
        alpha.y *= len;
        return alpha;
    }
#else
    /**
     * Get alpha vector (single precision, with scaling). For scaling
     * see setScale().
     *
     * @param theta Theta vector.
     * @returns Alpha vector.
     */
    __host__ __device__ Vector2D<float>
    getAlpha(const Vector2D<float> &theta) const {
        // printf("Scale factor plummer %f\n", m_scale);
        auto alpha = theta;
        alpha -= m_position;
        alpha /= (alpha.lengthSq() + m_angularwidth2);
        alpha *= m_4GM_f;
        return alpha;
    }
#endif
};
