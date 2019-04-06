#pragma once

#include "util/constants.h"
#include "util/vector2d.h"
#include <cstdio>

/**
 * Plummer lens class.
 */
class Plummer {
  public:
#ifdef __CUDACC__
    float4 m_data;
#else
    float m_angularwidth2;
    float m_4GM_f;
    Vector2D<float> m_position;
    char __padding[12];
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
        m_data.x = angularwidth * angularwidth;
        m_data.y = _4GM;
        m_data.z = position.x;
        m_data.w = position.y;
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
#ifdef __CUDACC__
        m_data.y = m_4GM * scalar;
#else
        m_4GM_f = m_4GM * scalar;
#endif
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
        float2 alpha = theta;
        alpha.x -= m_data.z;
        alpha.y -= m_data.w;
        float len = (alpha.x * alpha.x) + (alpha.y * alpha.y);
        len += m_data.x;
        len = (1 / len) * m_data.y;
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
