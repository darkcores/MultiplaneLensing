#pragma once

#include <cstdio>
#include "util/vector2d.h"
// CUDA includes
#include <cuda_runtime_api.h>

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
    __device__ Plummer() {}
    __host__ __device__ Plummer(const double Dd, const double mass,
                                const double angularwidth);

    __host__ __device__ Vector2D<double> getAlpha(Vector2D<double> theta) const;
    __host__ __device__ Vector2D<double> getBeta(Vector2D<double> theta) const;

    __host__ __device__ Vector2D<float>
    getAlphaf(const Vector2D<float> &theta) const {
		// printf("Scale factor plummer %f\n", m_scale);
        auto alpha = theta / (theta.lengthSq() + m_angularwidth2_f);
        alpha *= m_4GM_f;
        return alpha;
    }

    __host__ __device__ Vector2D<float>
    getBetaf(const Vector2D<float> &theta) const {
        auto beta = theta / (theta.lengthSq() + m_angularwidth2_f);
        beta *= m_4GM_s_f;
        beta = theta - beta;
        return beta;
    }

    __host__ __device__ Vector2D<float> getBetaf(const Vector2D<float> &theta,
                                                 const float &Ds,
                                                 const float &Dds) const {
        auto beta = theta / (theta.lengthSq() + m_angularwidth2_f);
        beta *= (Dds / Ds) * m_4GM_f;
        beta = theta - beta;
        return beta;
    }

    __host__ __device__ void setDistance(const double Dd);
    __host__ __device__ void setSource(const double Ds, const double Dds);
    __host__ __device__ void setScale(const float scale);
    __host__ __device__ void setMass(const double m_mass);
};
