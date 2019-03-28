#include "gtest/gtest.h"

#include <composite.h>
#include <util/constants.h>
#include <util/cosmology.h>

#include <iostream>

CompositeLensBuilder createGrid(double Dd, int N, double width, double height,
                                double angularwidth, double mass, double Ds = 0,
                                double Dds = 0) {
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
    lensbuilder.setSource(Ds, Dds);
    return lensbuilder;
}

TEST(CompositeTests, TestAlphaf) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d = 0.4;
    auto Dd = cosm.angularDiameterDistance(z_d);
    auto lensbuilder = createGrid(Dd, 3, 15 * ANGLE_ARCSEC, 15 * ANGLE_ARCSEC,
                                  5 * ANGLE_ARCSEC, 1e13 * MASS_SOLAR);
    lensbuilder.setScale(60);
    auto lens = lensbuilder.getLens();
    Vector2D<float> point(1 * ANGLE_ARCSEC, 2 * ANGLE_ARCSEC);
    auto alpha = lens.getAlphaf(point);
    // EXPECT_EQ(alpha.x(), 2.01259882e-05);
    // EXPECT_EQ(alpha.y(), 3.91372304e-05);
    EXPECT_LT(abs(alpha.x() - 2.01259882e-05), 1e-10);
    EXPECT_LT(abs(alpha.y() - 3.91372304e-05), 1e-10);
    lens.destroy();
}

TEST(CompositeTests, TestBetaf) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d = 0.4;
    double z_s = 1.0;
    auto Dd = cosm.angularDiameterDistance(z_d);
    auto Ds = cosm.angularDiameterDistance(z_s);
    auto Dds = cosm.angularDiameterDistance(z_d, z_s);
    auto lensbuilder = createGrid(Dd, 3, 15 * ANGLE_ARCSEC, 15 * ANGLE_ARCSEC,
                                  5 * ANGLE_ARCSEC, 1e13 * MASS_SOLAR, Ds, Dds);
    lensbuilder.setScale(30);
    lensbuilder.setSource(Ds, Dds);
    auto lens = lensbuilder.getLens();
    Vector2D<float> point(1 * ANGLE_ARCSEC, 2 * ANGLE_ARCSEC);
    // for (int i = 0; i < 16777216; i++) {
    auto beta = lens.getBetaf(point);
    EXPECT_LT(abs(beta.x() + 5.82627194e-06), 1e-11);
    EXPECT_LT(abs(beta.y() + 1.10613056e-05), 1e-11);
    // }
    // EXPECT_EQ(beta.x(), -5.82627194e-06);
    // EXPECT_EQ(beta.y(), -1.10613056e-05);
    lens.destroy();
}

// Test for some API wierdness
TEST(CompositeTests, TestBeta11) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d = 1.5;
    double z_s = 2.0;
    auto Dd = cosm.angularDiameterDistance(z_d);
    auto Ds = cosm.angularDiameterDistance(z_s);
    auto Dds = cosm.angularDiameterDistance(z_d, z_s);
    CompositeLensBuilder lensbuilder;

    Plummer plum(Dd, 1e13 * MASS_SOLAR * 11, 30 * ANGLE_ARCSEC);
    lensbuilder.addLens(plum, Vector2D<float>(0, 0));

    lensbuilder.setScale(3600);
    lensbuilder.setSource(Ds, Dds);
    auto lens = lensbuilder.getLens();
    Vector2D<float> point(1 * ANGLE_ARCSEC, 1 * ANGLE_ARCSEC);
    auto beta = lens.getBetaf(point);

    // EXPECT_EQ(beta.x(), 4.41364469e-6);
    // EXPECT_EQ(beta.y(), 4.41364469e-6);
    EXPECT_LT(abs(beta.x() - 4.41364469e-6), 1e-12);
    EXPECT_LT(abs(beta.y() - 4.41364469e-6), 1e-12);
	
    lens.destroy();
}
