#include "gtest/gtest.h"

#include <context.h>
#include <util/constants.h>

#include <fstream>
#include <iostream>

TEST(APITests, SetupTest) {
    // Setup variables
    const double unit = ANGLE_ARCSEC;
    const Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    std::vector<float> l_z = {1.5, 2.3};
    std::vector<float> s_z = {2.0, 2.7};
    std::vector<std::vector<PlummerParams>> params;
	
    for (auto x : l_z) {
        std::vector<PlummerParams> p;
        for (int i = 0; i < 10; i++) {
            PlummerParams pp = {Vector2D<float>(2, 7), 30.0,
                                10e13 * MASS_SOLAR};
            p.push_back(pp);
        }
        params.push_back(p);
    }

    std::vector<std::vector<Vector2D<float>>> thetas(10);
    for (int j = 0; j < 10; j++) {
        for (int i = 0; i < 100; i++) {
            thetas[j].push_back(Vector2D<float>(i, i + j));
        }
    }

    // Setup context
    MultiPlaneContext ctx(unit, cosm);
    ctx.init(l_z, params, s_z);
    ctx.setThetas(thetas);
}

/*
TEST(APITests, CalcTest) {
    // Setup variables
    const double unit = ANGLE_ARCSEC;
    const Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    std::vector<float> l_z = {1.5, 2.3};
    std::vector<float> s_z = {2.0, 2.7};
    std::vector<std::vector<PlummerParams>> params;
    for (auto x : l_z) {
        std::vector<PlummerParams> p;
        for (int i = 0; i < 1; i++) {
            PlummerParams pp = {Vector2D<float>(0, 0), 30.0,
                                10e15 * MASS_SOLAR};
            p.push_back(pp);
        }
        params.push_back(p);
    }

    std::vector<std::vector<Vector2D<float>>> thetas(2);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 10; j++) {
            thetas[i].push_back(Vector2D<float>(j, j));
        }
    }

    // Setup context
    MultiPlaneContext ctx(unit, cosm);
    ctx.init(l_z, params, s_z);
    ctx.setThetas(thetas);

    std::vector<std::vector<float>> masses(2);
    for (int x = 0; x < 2; x++) {
        for (int i = 0; i < 1; i++) {
            masses[x].push_back(11 - i);
        }
    }
    ctx.calculatePositions(masses);

    auto beta0 = ctx.getSourcePositions(0);

    EXPECT_LT(fabs(beta0[0].x() - 0.000000), 1e-6);
    EXPECT_LT(fabs(beta0[0].y() - 0.000000), 1e-6);
    EXPECT_LT(fabs(beta0[1].x() - 0.910380), 1e-6);
    EXPECT_LT(fabs(beta0[1].y() - 0.910380), 1e-6);
    EXPECT_LT(fabs(beta0[2].x() - 1.821944), 1e-6);
    EXPECT_LT(fabs(beta0[2].y() - 1.821944), 1e-6);
    EXPECT_LT(fabs(beta0[3].x() - 2.735825), 1e-6);
    EXPECT_LT(fabs(beta0[3].y() - 2.735825), 1e-6);
    EXPECT_LT(fabs(beta0[4].x() - 3.653057), 1e-6);
    EXPECT_LT(fabs(beta0[4].y() - 3.653057), 1e-6);
    EXPECT_LT(fabs(beta0[5].x() - 4.574539), 1e-6);
    EXPECT_LT(fabs(beta0[5].y() - 4.574539), 1e-6);
    EXPECT_LT(fabs(beta0[6].x() - 5.501002), 1e-6);
    EXPECT_LT(fabs(beta0[6].y() - 5.501002), 1e-6);
    EXPECT_LT(fabs(beta0[7].x() - 6.433003), 1.5e-6);
    EXPECT_LT(fabs(beta0[7].y() - 6.433003), 1.5e-6);
    EXPECT_LT(fabs(beta0[8].x() - 7.370913), 1.5e-6);
    EXPECT_LT(fabs(beta0[8].y() - 7.370913), 1.5e-6);
    EXPECT_LT(fabs(beta0[9].x() - 8.314935), 1.5e-6);
    EXPECT_LT(fabs(beta0[9].y() - 8.314935), 1.5e-6);

    auto beta1 = ctx.getSourcePositions(1);

    EXPECT_LT(fabs(beta1[0].x() - 0.000000), 1.25e-6);
    EXPECT_LT(fabs(beta1[0].y() - 0.000000), 1.25e-6);
    EXPECT_LT(fabs(beta1[1].x() - 0.836022), 1.25e-6);
    EXPECT_LT(fabs(beta1[1].y() - 0.836022), 1.25e-6);
    EXPECT_LT(fabs(beta1[2].x() - 1.674018), 1.25e-6);
    EXPECT_LT(fabs(beta1[2].y() - 1.674018), 1.25e-6);
    EXPECT_LT(fabs(beta1[3].x() - 2.515884), 1.25e-6);
    EXPECT_LT(fabs(beta1[3].y() - 2.515884), 1.25e-6);
    EXPECT_LT(fabs(beta1[4].x() - 3.363376), 1.25e-6);
    EXPECT_LT(fabs(beta1[4].y() - 3.363376), 1.25e-6);
    EXPECT_LT(fabs(beta1[5].x() - 4.218047), 1.25e-6);
    EXPECT_LT(fabs(beta1[5].y() - 4.218047), 1.25e-6);
    EXPECT_LT(fabs(beta1[6].x() - 5.081204), 1.25e-6);
    EXPECT_LT(fabs(beta1[6].y() - 5.081204), 1.25e-6);
    EXPECT_LT(fabs(beta1[7].x() - 5.953879), 1.25e-6);
    EXPECT_LT(fabs(beta1[7].y() - 5.953879), 1.25e-6);
    EXPECT_LT(fabs(beta1[8].x() - 6.836822), 1.25e-6);
    EXPECT_LT(fabs(beta1[8].y() - 6.836822), 1.25e-6);
    EXPECT_LT(fabs(beta1[9].x() - 7.730499), 1.25e-6);
    EXPECT_LT(fabs(beta1[9].y() - 7.730499), 1.25e-6);
}
*/

TEST(APITests, TestGrale1) {
    const double unit = ANGLE_ARCSEC;
    const Cosmology cosm(0.7, 0.3, 0.0, 0.7);
    std::vector<float> l_z = {0.4, 1.0};
    std::vector<float> s_z = {0.8, 1.5, 2.5};
    std::vector<std::vector<PlummerParams>> params;

    std::ifstream lensdata(
        "/home/jorrit/Sync/Universiteit/BAPR/lenstest/tests/testdata.txt");

    if (!lensdata.is_open())
        throw(1);

    std::vector<PlummerParams> p;
    for (int i = 0; i < 9; i++) {
        double x, y;
        lensdata >> x >> y;
        x /= unit;
        y /= unit;
        PlummerParams pp = {Vector2D<float>(x, y), 5.0, 1e13 * MASS_SOLAR};
        p.push_back(pp);
    }
    params.push_back(p);
    std::vector<PlummerParams> p2;
    for (int i = 0; i < 9; i++) {
        double x, y;
        lensdata >> x >> y;
        x /= unit;
        y /= unit;
        PlummerParams pp = {Vector2D<float>(x, y), 6.0, 1e13 * MASS_SOLAR};
        p2.push_back(pp);
    }
    params.push_back(p2);

    lensdata.close();

    std::ifstream thetadata(
        "/home/jorrit/Sync/Universiteit/BAPR/lenstest/tests/test_thetas.txt");

    std::vector<std::vector<Vector2D<float>>> thetas;
    int numthetas = 0;
    int t = 0;
    do {
        thetadata >> numthetas;
        // std::cout << "Next " << numthetas << " thetas" << std::endl;
        if (numthetas == 0)
            break;
        thetas.push_back(std::vector<Vector2D<float>>(numthetas));
        for (int i = 0; i < numthetas; i++) {
            double x, y;
            thetadata >> x >> y;
            x /= unit;
            y /= unit;
            thetas[t][i] = Vector2D<float>(x, y);
        }
        t++;
    } while (numthetas != 0);

    thetadata.close();

    // Setup context
    MultiPlaneContext ctx(unit, cosm);
    ctx.init(l_z, params, s_z);
    ctx.setThetas(thetas);

    std::vector<std::vector<float>> masses{{0, 1, 0, 1, 0, 2, 0, 0.5, 0},
                                           {1, 1, 0, 0.5, 1.5, 1, 0, 0, 0.5}};
    ctx.calculatePositions(masses);

    std::ifstream betadata(
        "/home/jorrit/Sync/Universiteit/BAPR/lenstest/tests/test_betas.txt");

    int numbetas = 0;
    betadata >> numbetas;
    t = 0;
    while (numbetas != 0) {
        const auto betavecs = ctx.getSourcePositions(t);
        // std::cout << "Next " << numbetas << " betas" << std::endl;
        ASSERT_EQ(numbetas, betavecs.size());
        for (int i = 0; i < numbetas; i++) {
            double x, y;
            betadata >> x >> y;
            // std::cout << "Also: " << thetas[i][0].x() << std::endl;
            // std::cout << "Also: " << x << std::endl;
            x /= unit;
            y /= unit;

            // printf("\r> %d i: %d, x: %f <> %f <> %f\t\t", t, i,
            //       betavecs[i].x(), x, thetas[t][i].x());
            ASSERT_LT(fabs(betavecs[i].x() - x), 1e-5);
            ASSERT_LT(fabs(betavecs[i].y() - y), 1e-5);
        }
        betadata >> numbetas;
        t++;
    }
	// printf("\n");

    betadata.close();
}