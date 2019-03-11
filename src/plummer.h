#pragma once

#include <cstdio>
#include "util/vector2d.h"
// CUDA includes
#include <cuda_runtime_api.h>

/**
 * Plummer lens class.
 */
class Plummer {
  private:
    double m_Dd, m_Ds, m_Dds;
    double m_mass;
    double m_angularwidth, m_angularwidth2;
    double m_4GM, m_4GM_s;
    // Vector2D<double> m_angularpos;
    float m_4GM_f, m_4GM_s_f, m_angularwidth2_f;
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
	 * Get alpha vector (double precision).
	 *
	 * @param theta Theta vector.
	 * @returns Alpha vector.
	 */
    __host__ __device__ Vector2D<double> getAlpha(Vector2D<double> theta) const;
	/**
	 * Get beta vector (double precision). You need to set the source
	 * first with setSource().
	 *
	 * @param theta Theta vector.
	 * @returns Beta vector.
	 */
    __host__ __device__ Vector2D<double> getBeta(Vector2D<double> theta) const;

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
        auto alpha = theta / (theta.lengthSq() + m_angularwidth2_f);
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
        auto beta = theta / (theta.lengthSq() + m_angularwidth2_f);
        beta *= m_4GM_s_f;
        beta = theta - beta;
        return beta;
    }

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
        auto beta = theta / (theta.lengthSq() + m_angularwidth2_f);
        beta *= (Dds / Ds) * m_4GM_f;
        beta = theta - beta;
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
};
