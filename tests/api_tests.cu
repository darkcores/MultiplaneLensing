#include "gtest/gtest.h"

#include <multiplane_cu.h>
#include <util/constants.h>

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

    std::vector<Vector2D<float>> thetas;
    for (int i = 0; i < 100; i++) {
        thetas.push_back(Vector2D<float>(i, i));
    }

    // Setup context
    MultiPlaneContext ctx(unit, cosm);
    ctx.init(l_z, params, s_z);
    ctx.setThetas(thetas);
}

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

    std::vector<Vector2D<float>> thetas;
    for (int i = 0; i < 10; i++) {
        thetas.push_back(Vector2D<float>(i, i));
    }

    // Setup context
    MultiPlaneContext ctx(unit, cosm);
    ctx.init(l_z, params, s_z);
    ctx.setThetas(thetas);

    std::vector<std::vector<double>> masses(2);
    for (int x = 0; x < 2; x++) {
        for (int i = 0; i < 1; i++) {
            masses[x].push_back(1e13 * MASS_SOLAR * (11 - i));
        }
    }
    ctx.calculatePositions(masses);

    auto beta0 = ctx.getSourcePositions(0);
    for (auto &x : beta0) {
        std::cout << "[ " << x.x() << " ; " << x.y() << " ]" << std::endl;
    }
    auto beta1 = ctx.getSourcePositions(1);
    for (auto &x : beta1) {
        std::cout << "[ " << x.x() << " ; " << x.y() << " ]" << std::endl;
    }
}
