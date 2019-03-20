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
    plum.setSource(Ds, Dds);
    auto trace = plum.getBeta(point);
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
    auto alpha = plum.getAlpha(point);
    // EXPECT_EQ(alpha.x(), 1.26193981e-41);
    // EXPECT_EQ(alpha.y(), 6.30969907e-42);
    EXPECT_LT(abs(alpha.x() - 1.26193981e-41), 1e-49);
    EXPECT_LT(abs(alpha.y() - 6.30969907e-42), 1e-50);
}

TEST(PlummerTests, TestAlphaF) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d = 1.5;
    auto Dd = cosm.angularDiameterDistance(z_d);
    Plummer plum(Dd, 1e5 * MASS_SOLAR, 30 * ANGLE_ARCSEC);
    plum.setScale(3600);
    Vector2D<float> point(1 * ANGLE_ARCSEC, 2 * ANGLE_ARCSEC);
    Vector2D<double> refpoint(1 * ANGLE_ARCSEC, 2 * ANGLE_ARCSEC);
    point *= (3600);
    auto alpha = plum.getAlphaf(point);
    alpha /= (3600);
    auto refalpha = plum.getAlpha(refpoint);
    // EXPECT_EQ(alpha.x(), refalpha.x());
    // EXPECT_EQ(alpha.y(), refalpha.y());
    EXPECT_LT(abs(alpha.x() - refalpha.x()), 1e-20);
    EXPECT_LT(abs(alpha.y() - refalpha.y()), 1e-20);
}

TEST(PlummerTests, TestBetaF) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d = 1.5;
    double z_s = 2.5;
    auto Dd = cosm.angularDiameterDistance(z_d);
    auto Ds = cosm.angularDiameterDistance(z_s);
    auto Dds = cosm.angularDiameterDistance(z_d, z_s);
    Plummer plum(Dd, 100000000000000, 30 * ANGLE_ARCSEC);
    plum.setSource(Ds, Dds);
    plum.setScale(60 * 60);
    Vector2D<float> point(1 * ANGLE_ARCMIN, 2 * ANGLE_ARCMIN);
    point *= (60 * 60);
    Vector2D<double> rpoint(1 * ANGLE_ARCMIN, 2 * ANGLE_ARCMIN);
    auto trace = plum.getBetaf(point);
    trace /= (60 * 60);
    // auto rtrace = plum.getBeta(rpoint);
    // EXPECT_EQ(trace.x(), rtrace.x());
    // EXPECT_EQ(trace.y(), rtrace.y());
    EXPECT_LT(abs(trace.x() - 0.00029089), 1e-8);
    EXPECT_LT(abs(trace.y() - 0.00058178), 1e-8);
}

// To look for a wierd result in the API
TEST(PlummerTests, TestBeta11) {
    Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    double z_d = 1.5;
    double z_s = 2.0;
    auto Dd = cosm.angularDiameterDistance(z_d);
    auto Ds = cosm.angularDiameterDistance(z_s);
    auto Dds = cosm.angularDiameterDistance(z_d, z_s);
    Plummer plum(Dd, 1e13 * MASS_SOLAR * 11, 30 * ANGLE_ARCSEC);
    plum.setSource(Ds, Dds);
    plum.setScale(60 * 60);
    Vector2D<float> point(1 * ANGLE_ARCSEC, 1 * ANGLE_ARCSEC);
    point *= (60 * 60);
    auto trace = plum.getBetaf(point);
    trace /= (60 * 60);
    // EXPECT_EQ(trace.x(), 4.41364469e-6);
    // EXPECT_EQ(trace.y(), 4.41364469e-6);
    EXPECT_LT(abs(trace.x() - 4.41364469e-6), 1e-12);
    EXPECT_LT(abs(trace.y() - 4.41364469e-6), 1e-12);
}
