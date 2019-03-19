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
