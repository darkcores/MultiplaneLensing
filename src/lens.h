#pragma once

#include "util/vector2d.h"
// CUDA includes
#include <cuda_runtime_api.h>

class Lens {
  private:
  public:
    __host__ __device__ virtual Vector2D<double>
    traceTheta(double Ds, double Dds, Vector2D<double> theta) const = 0;
    __host__ __device__ virtual Vector2D<double>
    getAlphaVector(Vector2D<double> theta) const = 0;
    __host__ __device__ virtual void
    getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy,
                              double &axy) const = 0;
    __host__ __device__ virtual double
    getSurfaceMassDensity(Vector2D<double> theta) const = 0;
};
