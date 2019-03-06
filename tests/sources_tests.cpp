#include "gtest/gtest.h"

#include <composite.h>
#include <multiplane.h>
#include <util/cosmology.h>
#include <util/constants.h>

#include <vector>
#include <iostream>

CompositeLensBuilder createGrid(double Dd, int N, double width, double height,
                                double angularwidth, double mass) {
    CompositeLensBuilder lensbuilder;
    double xstart = -width / 2;
    double xend = width / 2;
    double xstep = width / (N - 1);
    double ystart = -height / 2;
    double yend = height / 2;
    double ystep = height / (N - 1);
    // This works for this test, but might not always work TODO
    for (double x = xstart; x <= xend; x += xstep) {
        for (double y = ystart; y <= yend; y += ystep) {
            Plummer plum(Dd, mass, angularwidth);
            lensbuilder.addLens(plum, Vector2D<float>(x, y));
        }
    }
    return lensbuilder;
}


TEST(MultiplaneTests, TestBetaf) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d1 = 0.4;
	double z_d2 = 0.8;
    double z_s = 1.2;
    auto Dd1 = cosm.angularDiameterDistance(z_d1);
    auto Dd2 = cosm.angularDiameterDistance(z_d2);

    auto lensbuilder = createGrid(Dd1, 3, 15 * ANGLE_ARCSEC, 15 * ANGLE_ARCSEC,
                                  5 * ANGLE_ARCSEC, 1e13 * MASS_SOLAR);
	lensbuilder.setRedshift(z_d1);
    // auto lens1 = lensbuilder.getLens();
    auto lensbuilder2 = createGrid(Dd2, 3, 15 * ANGLE_ARCSEC, 15 * ANGLE_ARCSEC,
                                  5 * ANGLE_ARCSEC, 1e13 * MASS_SOLAR);
	lensbuilder2.setRedshift(z_d2);
	// auto lens2 = lensbuilder.getLens();

	MultiplaneBuilder planebuilder(cosm);
	planebuilder.addPlane(&lensbuilder);
	planebuilder.addPlane(&lensbuilder2);

	SourcePlaneBuilder sourcebuilder(z_s);
	sourcebuilder.addPoint(Vector2D<float>(1 * ANGLE_ARCSEC, 1 * ANGLE_ARCSEC), 1 * ANGLE_ARCSEC);
	auto sourceplane = sourcebuilder.getPlane();
	planebuilder.addSourcePlane(sourceplane);

	auto multiplane = planebuilder.getMultiPlane();

	std::cout << "Setup done" << std::endl;
	
    Vector2D<float> point(1 * ANGLE_ARCSEC, 1 * ANGLE_ARCSEC);
	auto pixel = multiplane.traceTheta(point);
	ASSERT_EQ(pixel, 0);
    // for (int i = 0; i < 16777216; i++) {
    // auto beta = lens.getBetaf(point, Ds, Dds);
    // EXPECT_LT(abs(beta.x() + 5.82627194e-06), 1e-11);
    // EXPECT_LT(abs(beta.y() + 1.10613056e-05), 1e-11);
    // }
    // EXPECT_EQ(beta.x(), -5.82627194e-06);
    // EXPECT_EQ(beta.y(), -1.10613056e-05);
}
