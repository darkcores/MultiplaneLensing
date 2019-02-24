#pragma once

#include "lens.h"
// CUDA includes
#include <cuda_runtime_api.h>

class Plummer : public Lens {
  private:
    double m_Dd;
    double m_mass;
    double m_angularwidth, m_angularwidth2;
    // Vector2D<double> m_angularpos;

  public:
    __host__ __device__ Plummer(const double Dd, const double mass, const double angularwidth/*,
																		   const Vector2D<double> &angularposition*/);

    __host__ __device__ virtual Vector2D<double>
    traceTheta(double Ds, double Dds, Vector2D<double> theta) const;
    __host__ __device__ virtual Vector2D<double>
    getAlphaVector(Vector2D<double> theta) const;
    __host__ __device__ virtual void
    getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy,
                              double &axy) const;
    __host__ __device__ virtual double
    getSurfaceMassDensity(Vector2D<double> theta) const;

  private:
    __host__ __device__ double getMassInside(double thetaLength) const;
    __host__ __device__ double getLensDistance() const { return m_Dd; }
    __host__ __device__ double
    getProfileSurfaceMassDensity(double thetaLength) const;
};
