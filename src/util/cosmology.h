#pragma once

// CUDA includes
#ifdef __CUDACC__
#include <cuda_runtime_api.h>
#else
#ifndef __host__
#define __host__
#define __device__
#endif
#endif

class Cosmology {
  private:
    double m_h, m_w;
    double m_Wm, m_Wr, m_Wv, m_Wk;

  public:
	/**
	 * Empty constructor, don't know why. Might delete later.
	 */
	__host__ __device__ Cosmology() {}
	/**
	 * Cosmology calculations class.
	 *
	 * @param h Hubble constant.
	 * @param omega_m Matter density parameter.
	 * @param omega_r Radiation density parameter.
	 * @param omega_v Vacuum density parameter.
	 * @param w Equation of state of the vacuum energy.
	 */
    __host__ __device__ Cosmology(const double h, const double omega_m,
                                  const double omega_r, const double omega_v,
                                  const double w = -1.0);

	/**
	 * Calculate angular diameter distance for redshift(s).
	 * TODO: Use GSL for integration
	 *
	 * @param z1 First redshift.
	 * @param z2 Second redshift.
	 * @param num
	 * @returns Angular diameter distance.
	 */
    __host__ __device__ double
    angularDiameterDistance(double z1, double z2 = -1.0, int num = 10000) const;

	/*
    __host__ __device__ double
    redshiftForAngularDiameterDistance(double Dtarget,
                                       double zref = 0,
                                       double zmax = 20) const;
	*/
};
