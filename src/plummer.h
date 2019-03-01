#pragma once

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
    // Scaling doesn't work here I think TODO
    __host__ __device__ Vector2D<float> getAlphaf(const Vector2D<float> &theta) const;
    __host__ __device__ Vector2D<double> getBeta(Vector2D<double> theta) const;
    __host__ __device__ Vector2D<float> getBetaf(const Vector2D<float> &theta) const;

    __host__ __device__ void setDistance(const double Dd);
    __host__ __device__ void setSource(const double Ds, const double Dds);
    __host__ __device__ void setScale(const float scale = 60);
    __host__ __device__ void setMass(const double m_mass);
};

/*
class Plummerf {
  private:
    float m_Dd, m_D;
	float m_mass, m_angularwidth;
    float m_4GM_f, m_angularwidth2_f;

  public:
    __host__ __device__ Plummerf(const float Dd, const float mass,
                                 const float angularwidth);

    __host__ __device__ Vector2D<float> getAlphaf(Vector2D<float> theta) const;
    __host__ __device__ Vector2D<float> getBetaf(Vector2D<float> theta) const;

    __host__ __device__ void setMass(const float m_mass);
    __host__ __device__ void setSource(const float Ds, const float Dds);
    __host__ __device__ void setScale(const float scale = 60);
};
*/
