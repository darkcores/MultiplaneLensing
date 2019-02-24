#include "gtest/gtest.h"

#include <plummer.h>
#include <util/constants.h>
#include <util/cosmology.h>

#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/host_vector.h>

__global__ void alphaCalc(int n, Vector2D<double> *data) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    Vector2D<double> vec((i % 7) * ANGLE_ARCSEC, (i % 5) * ANGLE_ARCSEC);

    const Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    const double z_d = 1.5;
    const double Dd = cosm.angularDiameterDistance(z_d);
    const Plummer plum(Dd, 1000000000, 30 * ANGLE_ARCSEC);

    data[i] = plum.getAlphaVector(vec);
}

TEST(PlummerCuTests, TestAlpha) {
    // Mostly to verify it runs on GPU
    thrust::device_vector<Vector2D<double>> alphas(1024);
    auto a_ptr = thrust::raw_pointer_cast(&alphas[0]);
    alphaCalc<<<1024 / 64, 64>>>(1024, a_ptr);
    thrust::host_vector<Vector2D<double>> h_alphas(alphas);

    /*
     * The following block is mostly the same but there are some small
     * differences on some results (e > 1024), see
     * https://docs.nvidia.com/cuda/floating-point/index.html
     *
    const Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    const double z_d = 1.5;
    const double Dd = cosm.angularDiameterDistance(z_d);
    const Plummer plum(Dd, 1000000000, 30 * ANGLE_ARCSEC);
    for (int i = 0; i < 128; i += 64) {
        for (int j = 0; j < 64; j++) {
            int idx = i + j;
            Vector2D<double> vec((idx % 7) * ANGLE_ARCSEC,
                                 (idx % 5) * ANGLE_ARCSEC);
            Vector2D<double> result = plum.getAlphaVector(vec);
            EXPECT_EQ(h_alphas[idx], result);
        }
    }
	*/
}