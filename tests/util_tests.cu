#include "gtest/gtest.h"

#include <util/constants.h>
#include <util/cosmology.h>
#include <util/vector2d.h>
#include <util/error.h>

#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/host_vector.h>

__global__ void scalegpu(int n, Vector2D<double> *data, double scalar) {
    int i = threadIdx.x;

    if (i < n) {
        data[i].setX(i);
        data[i].setY(i);
        data[i] *= scalar;
    }
}

__global__ void cosmologycalc(int n, double *redshifts, double *angdists) {
    int i = threadIdx.x;
    if (i < n) {
        const Cosmology cosm(0.7, 0.3, 0.0, 0.7);
        angdists[i] = cosm.angularDiameterDistance(0.0, redshifts[i]);
    }
}

template <typename T> struct mul {
    __device__ T operator()(const T &x, double y) const { return x * y; }
};

TEST(UtilCuTests, ErrTests) {
	try {
		gpuErrchk((cudaError_t)-1);
	} catch (int e) {
		ASSERT_EQ(e, -1);
	}
}

TEST(UtilCuTests, ScalingDouble) {
    thrust::host_vector<Vector2D<double>> H(4);
    thrust::device_vector<Vector2D<double>> D = H;
    auto raw = thrust::raw_pointer_cast(&D[0]);
    scalegpu<<<1, 32>>>(4, raw, 4.0);
    mul<Vector2D<double>> op;
    thrust::transform(D.begin(), D.end(), thrust::make_constant_iterator(0.25),
                      D.begin(), op);
    thrust::copy(D.begin(), D.end(), H.begin());
    ASSERT_EQ(H[0].x(), 0);
    ASSERT_EQ(H[1].x(), 1);
}

TEST(UtilCuTests, CosmologyTests) {
    thrust::device_vector<double> redshifts(128);
    thrust::device_vector<double> angdists(128);
    thrust::sequence(redshifts.begin(), redshifts.end());
    auto r_ptr = thrust::raw_pointer_cast(&redshifts[0]);
    auto a_ptr = thrust::raw_pointer_cast(&angdists[0]);
    cosmologycalc<<<1, 128>>>(128, r_ptr, a_ptr);
    const double x1 = angdists[1] / DIST_MPC;
    ASSERT_EQ((int)x1, 1651);
    const double x127 = angdists[127] / DIST_MPC;
    ASSERT_EQ((int)x127, 99);
}
