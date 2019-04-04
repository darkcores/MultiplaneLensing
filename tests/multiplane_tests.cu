#include "gtest/gtest.h"

#include <multiplane.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

CompositeLensBuilder createGrid(float z, double Dd, int N, double width,
                                double height, double angularwidth,
                                double mass) {
    CompositeLensBuilder lensbuilder(z);
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
        }
    }
    return lensbuilder;
}

TEST(MultiplaneTests, TestTrace) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d1 = 0.4;
    double z_d2 = 1.0;
    auto Dd1 = cosm.angularDiameterDistance(z_d1);
    auto Dd2 = cosm.angularDiameterDistance(z_d2);

    MultiplaneBuilder builder(cosm);
    auto lensbuilder = createGrid(z_d1, Dd1, 3, 15, 15, 5, 1e13 * MASS_SOLAR);
    builder.addPlane(lensbuilder);
    auto lensbuilder2 = createGrid(z_d2, Dd2, 3, 15, 15, 5, 1e13 * MASS_SOLAR);
    builder.addPlane(lensbuilder2);

    std::vector<float> srcs{0.2, 0.8, 1.5, 2.5};
    builder.setRedshifts(srcs);

    auto mp = builder.getCuMultiPlane();

    float2 point{.x = 2, .y = 1};
	thrust::host_vector<float2> thetas;
	thetas.push_back(point);
	thrust::host_vector<float2> betas(1);

	thrust::device_vector<float2> d_thetas(thetas);
	thrust::device_vector<float2> d_betas(1);

	float2 *thetasptr = thrust::raw_pointer_cast(&d_thetas[0]);
	float2 *betasptr = thrust::raw_pointer_cast(&d_betas[0]);

	mp.traceThetas(thetasptr, betasptr, 1, 0);
	betas = d_betas;
	EXPECT_EQ(betas[0].x, 2.0);
	EXPECT_EQ(betas[0].y, 1.0);
	// printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

	mp.traceThetas(thetasptr, betasptr, 1, 1);
	betas = d_betas;
	EXPECT_LT(fabs(betas[0].x + 1.579611), 1e-5);
	EXPECT_LT(fabs(betas[0].y + 0.840784), 1e-5);
	// printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

	mp.traceThetas(thetasptr, betasptr, 1, 2);
	betas = d_betas;
	EXPECT_LT(fabs(betas[0].x + 1.736920), 1e-5);
	EXPECT_LT(fabs(betas[0].y + 0.876505), 1e-5);
	// printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

	mp.traceThetas(thetasptr, betasptr, 1, 3);
	betas = d_betas;
	EXPECT_LT(fabs(betas[0].x + 1.310565), 1e-5);
	EXPECT_LT(fabs(betas[0].y + 0.621893), 1e-5);
	// printf("Beta: %f; %f\n", betas[0].x(), betas[0].y());

    mp.destroy();
}