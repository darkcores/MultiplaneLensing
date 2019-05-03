#pragma once

#include "util/constants.h"
#include "util/vector2d.h"
#include <cstdio>

/**
 * Plummer lens class.
 */
class Plummer {
  public:
    float m_angularwidth2;
    float m_4GM_f;

#ifdef __CUDACC__
    float2 m_position;
#else
    Vector2D<float> m_position;
    char __padding[4];
#endif
    float m_4GM;

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
            const double scale, const float2 position) {
        float _4GM =
            (4 * CONST_G * mass) / (SPEED_C * SPEED_C * Dd) * (scale * scale);
        m_angularwidth2 = angularwidth * angularwidth;
        m_4GM_f = _4GM;
        m_position.x = position.x;
        m_position.y = position.y;
        m_4GM = _4GM;
        // printf("\x1B[34m4GM: %f (m %lf)\x1B[0m \n", m_4GM_f, mass);
    }

    __host__ __device__ Plummer() {}
#else
    Plummer(const double Dd, const double mass, const double angularwidth,
            const double scale, const Vector2D<float> position)
        : m_angularwidth2(angularwidth * angularwidth), m_position(position) {
        float _4GM =
            (4 * CONST_G * mass) / (SPEED_C * SPEED_C * Dd) * (scale * scale);
        m_4GM_f = _4GM;
        m_4GM = _4GM;
        // printf("\x1B[35m4GM: %f (m %lf)\x1B[0m \n", m_4GM_f, mass);
    }
#endif

    __host__ __device__ void update(const float scalar) {
        m_4GM_f = m_4GM * scalar;
        // printf("\x1B[31m4GM: %f (* %f)\x1B[0m \n", m_4GM_f, scalar);
    }

#ifdef __CUDACC__
    __device__ float4 f4() const {
        return float4{m_angularwidth2, m_4GM_f, m_position.x, m_position.y};
    }
    /**
     * Get alpha vector (single precision, with scaling). For scaling
     * see setScale().
     *
     * @param theta Theta vector.
     * @returns Alpha vector.
     */
    __host__ __device__ float2 getAlpha(const float2 &theta) const {
        // printf("Scale factor plummer %f\n", m_scale);
        float2 alpha = theta;
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
    inline __host__ __device__ Vector2D<float>
    getAlpha(const Vector2D<float> &theta) const {
        // printf("Scale factor plummer %f\n", m_scale);
        auto alpha = theta;
        alpha -= m_position;
		float len = alpha.lengthSq() + m_angularwidth2;
		len = (1 / len) * m_4GM_f;
        alpha *= len;
        return alpha;
    }
#endif
};
