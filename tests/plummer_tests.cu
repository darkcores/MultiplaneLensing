#include "gtest/gtest.h"

#include <plummer.h>
#include <util/constants.h>
#include <util/cosmology.h>

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

__global__ void alphaCalc(int n, float2 *thetas, float2 *betas, Plummer plum) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        // printf("Thread: %i - Theta: %f\n", i, thetas[i].x);
		betas[i] = plum.getAlpha(thetas[i]);
	}
}

TEST(PlummerCuTests, TestAlpha) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d = 1.5;
    double z_s = 2.0;
    auto Dd = cosm.angularDiameterDistance(z_d);
    auto Ds = cosm.angularDiameterDistance(z_s);
    auto Dds = cosm.angularDiameterDistance(z_d, z_s);
    float2 pos{.x = 0, .y = 0};
    thrust::host_vector<float2> thetas;
    Plummer plum(Dd, 1e13 * MASS_SOLAR, 30, 1 / ANGLE_ARCSEC, pos);
    for (int i = 0; i < 10; i++) {
        float2 point{.x = (float)i, .y = (float)i};
        thetas.push_back(point);
    }
    thrust::device_vector<float2> d_thetas(thetas), d_betas(10);

    float2 *thetaptr = thrust::raw_pointer_cast(&d_thetas[0]);
    float2 *betaptr = thrust::raw_pointer_cast(&d_betas[0]);
    alphaCalc<<<1, 32>>>(10, thetaptr, betaptr, plum);

    thrust::host_vector<float2> betas(d_betas);
    for (int i = 0; i < 10; i++) {
        betas[i].x = thetas[i].x - (betas[i].x * (Dds / Ds));
        betas[i].y = thetas[i].y - (betas[i].y * (Dds / Ds));
    }

    EXPECT_LT(fabs(betas[0].x - 0.0000000000), 1.25e-6);
    EXPECT_LT(fabs(betas[0].y - 0.0000000000), 1.25e-6);
    EXPECT_LT(fabs(betas[1].x - 0.9918526879), 1.25e-6);
    EXPECT_LT(fabs(betas[1].y - 0.9918526879), 1.25e-6);
    EXPECT_LT(fabs(betas[2].x - 1.9838130496), 1.25e-6);
    EXPECT_LT(fabs(betas[2].y - 1.9838130496), 1.25e-6);
    EXPECT_LT(fabs(betas[3].x - 2.9759840670), 1.25e-6);
    EXPECT_LT(fabs(betas[3].y - 2.9759840670), 1.25e-6);
    EXPECT_LT(fabs(betas[4].x - 3.9684597618), 1.25e-6);
    EXPECT_LT(fabs(betas[4].y - 3.9684597618), 1.25e-6);
    EXPECT_LT(fabs(betas[5].x - 4.9613217080), 1.25e-6);
    EXPECT_LT(fabs(betas[5].y - 4.9613217080), 1.25e-6);
    EXPECT_LT(fabs(betas[6].x - 5.9546365711), 1.25e-6);
    EXPECT_LT(fabs(betas[6].y - 5.9546365711), 1.25e-6);
    EXPECT_LT(fabs(betas[7].x - 6.9484547811), 1.25e-6);
    EXPECT_LT(fabs(betas[7].y - 6.9484547811), 1.25e-6);
    EXPECT_LT(fabs(betas[8].x - 7.9428103075), 1.25e-6);
    EXPECT_LT(fabs(betas[8].y - 7.9428103075), 1.25e-6);
    EXPECT_LT(fabs(betas[9].x - 8.9377213942), 1.25e-6);
    EXPECT_LT(fabs(betas[9].y - 8.9377213942), 1.25e-6);
}