#include "gtest/gtest.h"

#include <plummer.h>
#include <util/constants.h>
#include <util/cosmology.h>

TEST(PlummerTests, TestTraceTheta) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);

    // Test redshift values
    double z_d = 1.5;
    double z_s = 2.5;

    auto Dd = cosm.angularDiameterDistance(z_d);
	// EXPECT_EQ(Dd, 5.386180745148062e+25);
	EXPECT_LT(fabs(Dd - 5.386180745148062e+25), 1e+17);
    auto Ds = cosm.angularDiameterDistance(z_s);
    auto Dds = cosm.angularDiameterDistance(z_d, z_s);
	
	// Keep in mind the numbers we compare are a lot more limited in
	// precision.
    Plummer plum(Dd, 10000, 30 * ANGLE_ARCSEC);
    Vector2D<double> point(1 * ANGLE_ARCSEC, 0.5 * ANGLE_ARCSEC);
    auto trace = plum.traceTheta(Ds, Dds, point);
    // EXPECT_EQ(trace.x(), 4.84813681e-06);
    // EXPECT_EQ(trace.y(), 2.42406841e-06);
    EXPECT_LT(fabs(trace.x() - 4.84814e-06), 1e-11);
    EXPECT_LT(fabs(trace.y() - 2.42406e-06), 1e-11);
}

TEST(PlummerTests, TestAlpha) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);

    // Test redshift values
    double z_d = 1.5;

    auto Dd = cosm.angularDiameterDistance(z_d);
	EXPECT_LT(fabs(Dd - 5.386180745148062e+25), 1e+17);
    Plummer plum(Dd, 1000000000, 30 * ANGLE_ARCSEC);
    Vector2D<double> point(1 * ANGLE_ARCSEC, 0.5 * ANGLE_ARCSEC);
    auto alpha = plum.getAlphaVector(point);
    // EXPECT_EQ(alpha.x(), 1.26193981e-41);
    // EXPECT_EQ(alpha.y(), 6.30969907e-42);
    EXPECT_LT(abs(alpha.x() - 1.26193981e-41), 1e-49);
    EXPECT_LT(abs(alpha.y() - 6.30969907e-42), 1e-50);
}
