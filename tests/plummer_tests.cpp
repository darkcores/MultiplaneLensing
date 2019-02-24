#include "gtest/gtest.h"

#include <plummer.h>
#include <util/constants.h>
#include <util/cosmology.h>

TEST(PlummerTests, AlphaTheta) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);

    // Test redshift values
    double z_d = 1.5;
    double z_s = 2.5;

    auto Dd = cosm.angularDiameterDistance(z_d);
    auto Ds = cosm.angularDiameterDistance(z_s);
    auto Dds = cosm.angularDiameterDistance(z_d, z_s);

	Plummer plum(Dd, 10000, 30*ANGLE_ARCSEC);
	Vector2D<double> point(1*ANGLE_ARCSEC, 0.5*ANGLE_ARCSEC);
	auto trace = plum.traceTheta(Ds, Dds, point);
	ASSERT_LT(abs(trace.x() - 4.84814e-06), 0.000000001);
	ASSERT_LT(abs(trace.y() - 2.42406e-06), 0.000000001);
	auto alpha = plum.getAlphaVector(point);
	// ASSERT_EQ(alpha.x(), 1.74950205e-46);
	ASSERT_LT(abs(alpha.x() - 1.74950205e-46), 1e-45);
	ASSERT_LT(abs(alpha.y() - 8.74751025e-47), 1e-46);
}
