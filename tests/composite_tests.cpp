#include "gtest/gtest.h"

#include <composite.h>
#include <util/constants.h>
#include <util/cosmology.h>

#include <iostream>

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
                         Vector2D<float>(x, y));
            lensbuilder.addLens(plum);
        }
    }
    return lensbuilder;
}

TEST(CompositeTests, TestAlpha) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d = 0.4;
    auto Dd = cosm.angularDiameterDistance(z_d);
    auto lensbuilder = createGrid(Dd, 3, 15, 15, 5, 1e13 * MASS_SOLAR);
    auto lens = lensbuilder.getLens();
    Vector2D<float> point(1, 2);
    auto alpha = lens.getAlpha(point);
	alpha *= ANGLE_ARCSEC;
    // EXPECT_EQ(alpha.x(), 2.01259882e-05);
    // EXPECT_EQ(alpha.y(), 3.91372304e-05);
    EXPECT_LT(abs(alpha.x() - 2.01259882e-05), 1e-10);
    EXPECT_LT(abs(alpha.y() - 3.91372304e-05), 1e-10);
    lens.destroy();
}

TEST(CompositeTests, TestAlphaUpdate) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d = 0.4;
    double z_s = 0.8;
    auto Dd = cosm.angularDiameterDistance(z_d);
    auto Ds = cosm.angularDiameterDistance(z_s);
    auto Dds = cosm.angularDiameterDistance(z_d, z_s);
    auto lensbuilder = createGrid(Dd, 3, 15, 15, 5, 1e13 * MASS_SOLAR);
    auto lens = lensbuilder.getLens();

    // std::vector<float> factors{0, 1, 0, 1, 0, 2, 0, 0.5, 0};
    std::vector<float> factors{0, 1, 0, 1, 0, 0.5, 0, 2, 0};
	// std::reverse(factors.begin(), factors.end());
	lens.update(&factors[0]);
	
    Vector2D<float> point(2, 1);
    auto alpha = lens.getAlpha(point);
	auto trace = point - (alpha * (Dds / Ds));
    // EXPECT_EQ(trace.x(), 4.607982);
    // EXPECT_EQ(trace.y(), -1.687718);
    EXPECT_LT(abs(trace.x() - 4.607982), 1e-5);
    EXPECT_LT(abs(trace.y() + 1.687718), 1e-5);
    lens.destroy();}

