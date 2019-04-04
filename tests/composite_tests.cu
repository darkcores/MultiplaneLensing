#include "gtest/gtest.h"

#include <composite.h>
#include <thrust/device_vector.h>
#include <util/constants.h>
#include <util/cosmology.h>

CompositeLensBuilder createGrid(double Dd, int N, double width, double height,
                                double angularwidth, double mass) {
    CompositeLensBuilder lensbuilder(Dd); // Redshift not used in these tests.
    double xstart = -width / 2;
    double xend = width / 2;
    double xstep = width / (N - 1);
    double ystart = -height / 2;
    double yend = height / 2;
    double ystep = height / (N - 1);
    // This works for this test, but might not always work TODO
    for (double x = xstart; x <= xend; x += xstep) {
        for (double y = ystart; y <= yend; y += ystep) {
            Plummer plum(Dd, mass, angularwidth, 1 / ANGLE_ARCSEC,
                         float2{.x = (float)x, .y = (float)y});
            lensbuilder.addLens(plum);
			// printf("Added lens\n");
        }
    }
    return lensbuilder;
}

__global__ void alphaCompCalc(int n, float2 *thetas, float2 *alphas,
                              CompositeLens lens) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < n) {
        // printf("Thread: %i - Theta: %f\n", i, thetas[i].x);
        alphas[i] = lens.getAlpha(thetas[i]);
    }
}

TEST(CompositeTests, TestAlpha) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d = 0.4;
    auto Dd = cosm.angularDiameterDistance(z_d);
    auto lensbuilder = createGrid(Dd, 3, 15, 15, 5, 1e13 * MASS_SOLAR);
    auto lens = lensbuilder.getCuLens();

	thrust::host_vector<float2> thetas;
	thetas.push_back(float2{.x = 1.0, .y = 2.0});
    thrust::device_vector<float2> d_thetas(thetas), d_alphas(1);

    float2 *thetaptr = thrust::raw_pointer_cast(&d_thetas[0]);
    float2 *alphaptr = thrust::raw_pointer_cast(&d_alphas[0]);
    alphaCompCalc<<<1, 32>>>(1, thetaptr, alphaptr, lens);

	thrust::host_vector<float2> alphas(d_alphas);
    float2 alpha = alphas[0];
    alpha.x *= ANGLE_ARCSEC;
    alpha.y *= ANGLE_ARCSEC;
    // EXPECT_EQ(alpha.x, 2.01259882e-05);
    // EXPECT_EQ(alpha.y, 3.91372304e-05);
    EXPECT_LT(abs(alpha.x - 2.01259882e-05), 1e-10);
    EXPECT_LT(abs(alpha.y - 3.91372304e-05), 1e-10);
    lens.destroy();
}
