#pragma once

// CUDA includes
#include <cuda_runtime_api.h>

class Cosmology {
  private:
    double m_h, m_w;
    double m_Wm, m_Wr, m_Wv, m_Wk;

  public:
    __host__ __device__ Cosmology(const double h, const double omega_m,
                                  const double omega_r, const double omega_v,
                                  const double w = -1.0);

    __host__ __device__ double
    angularDiameterDistance(double z1, double z2 = -1.0, int num = 10000) const;

    __host__ __device__ double
    redshiftForAngularDiameterDistance(double Dtarget,
                                       double zref = 0,
                                       double zmax = 20) const;
};
