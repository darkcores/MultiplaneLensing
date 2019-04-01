#pragma once

#include "util/vector2d.h"
#include <cstdio>

class MiniPlummer {
  private:
    float m_angularwidth2_f, m_4GM_f;

  public:
    __host__ __device__ MiniPlummer() {}
    __host__ __device__ MiniPlummer(float GM, float angwidth) {
        m_angularwidth2_f = angwidth;
        m_4GM_f = GM;
    }

#ifdef __CUDACC__
    /**
     * Get alpha vector (single precision, with scaling). For scaling
     * see setScale().
     *
     * @param theta Theta vector.
     * @returns Alpha vector.
     */
    __host__ __device__ float2 getAlphaf(const float2 &theta) const {
        // printf("Scale factor plummer %f\n", m_scale);
        auto alpha = theta;
        float len = (theta.x * theta.x) + (theta.y * theta.y);
        len += m_angularwidth2_f;
        len = 1 / len;
        alpha.x *= len;
        alpha.y *= len;
        alpha.x *= m_4GM_f;
        alpha.y *= m_4GM_f;
        return alpha;
    }
#endif
};

/**
 * Plummer lens class.
 */
class Plummer {
  private:
    float m_Dd, m_Df;
    double m_mass;
    float m_angularwidth;
    // Vector2D<double> m_angularpos;
    float m_4GM_f, m_angularwidth2_f;
    float m_scale;

  public:
    /**
     * Create a new plummer lens, only on device.
     */
    __device__ Plummer() {}
    /**
     * Create a new plummer lens.
     *
     * @param Dd angular distance for this lens.
     * @param mass lens mass.
     * @param angularwidth angular width for this lens.
     */
    __host__ __device__ Plummer(const double Dd, const double mass,
                                const double angularwidth);

    /**
     * Get alpha vector (single precision, with scaling). For scaling
     * see setScale().
     *
     * @param theta Theta vector.
     * @returns Alpha vector.
     */
    __host__ __device__ Vector2D<float>
    getAlphaf(const Vector2D<float> &theta) const {
        // printf("Scale factor plummer %f\n", m_scale);
        auto alpha = theta;
        alpha /= (theta.lengthSq() + m_angularwidth2_f);
        alpha *= m_4GM_f;
        return alpha;
    }

    /**
     * Get beta vector (single precision, with scaling). For scaling
     * see setScale(). You need to set the source first with
     * setSource().
     *
     * @param theta Theta vector.
     * @returns Beta vector.
     */
    __host__ __device__ Vector2D<float>
    getBetaf(const Vector2D<float> &theta) const {
        auto beta = theta;
        beta /= (theta.lengthSq() + m_angularwidth2_f);
        beta *= m_4GM_f * m_Df;
        beta.rdiff(theta);
        return beta;
    }

#ifdef __CUDACC__
    /**
     * Get alpha vector (single precision, with scaling). For scaling
     * see setScale().
     *
     * @param theta Theta vector.
     * @returns Alpha vector.
     */
    __host__ __device__ float2 getAlphaf(const float2 &theta) const {
        // printf("Scale factor plummer %f\n", m_scale);
        auto alpha = theta;
        float len = (theta.x * theta.x) + (theta.y * theta.y);
        len += m_angularwidth2_f;
        len = 1 / len;
        alpha.x *= len;
        alpha.y *= len;
        alpha.x *= m_4GM_f;
        alpha.y *= m_4GM_f;
        return alpha;
    }
#endif

    /**
     * Get beta vector (single precision, with scaling). For scaling
     * see setScale().
     *
     * @param theta Theta vector.
     * @param Ds Source angular distance.
     * @param Dds lens<->source angular distance.
     * @returns Beta vector.
     */
    __host__ __device__ Vector2D<float> getBetaf(const Vector2D<float> &theta,
                                                 const float &Ds,
                                                 const float &Dds) const {
        auto beta = theta;
        beta /= (theta.lengthSq() + m_angularwidth2_f);
        beta *= (Dds / Ds) * m_4GM_f;
        beta.rdiff(theta);
        return beta;
    }

    /**
     * Set the angular distance for this lens.
     *
     * @param Dd angular distance.
     */
    __host__ __device__ void setDistance(const double Dd);

    /**
     * Set the source angular distance.
     *
     * @param Ds Source angular distance.
     * @param Dds lens<->source angular distance.
     */
    __host__ __device__ void setSource(const double Ds, const double Dds);

    /**
     * set single precision float scale factor.
     *
     * @param scale scale factor.
     */
    __host__ __device__ void setScale(const float scale);

    /**
     * Set lens mass.
     *
     * @param mass new lens mass.
     */
    __host__ __device__ void setMass(const double mass);

    __host__ __device__ MiniPlummer getMini() const {
        return MiniPlummer(m_4GM_f, m_angularwidth * m_angularwidth);
    }
};
