#include "gtest/gtest.h"

#include <composite.h>
#include <util/constants.h>
#include <util/cosmology.h>

#include <iostream>

CompositeLensBuilder createGrid(double Dd, int N, double width, double height,
                                double angularwidth, double mass, double Ds = 0,
                                double Dds = 0) {
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

TEST(CompositeTests, TestAlphaf) {
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
