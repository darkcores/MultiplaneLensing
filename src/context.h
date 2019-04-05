#pragma once

#include <vector>
// Cosmology calculations are needed anyways
#include "util/cosmology.h"
// Same for points/Vector2D
#include "util/vector2d.h"
// some template definitions
class Multiplane;
class CompositeLensBuilder;

/**
 * Initial parameters for each plummer lens.
 */
struct PlummerParams {
    Vector2D<float> position;
    float angularwidth;
    double mass;
};

/**
 * Context for multiplane plummer calculations.
 *
 * For sample use cases see api_tests.cu in tests/ or example.cpp in
 * example/ for a complete example program.
 */
class MultiPlaneContext {
  private:
    const double m_angularUnit;
    const Cosmology m_cosmology;
    Vector2D<float> *m_theta;
    size_t m_theta_len;
    std::vector<size_t> m_theta_count;
    Vector2D<float> *m_beta;
    Multiplane *m_multiplane;
    std::vector<std::vector<Vector2D<float>>> m_betas;

    CompositeLensBuilder buildLens(const float redshift,
                                   const std::vector<PlummerParams> &params);

  public:
    /**
     * Create a new context.
     *
     * @param angularUnit What unit is used in the input data,
     * eg. ANGLE_ARCSEC.
     * @param cosmology What parameters are used for cosmology.
     */
    MultiPlaneContext(const double angularUnit, const Cosmology cosmology);
    ~MultiPlaneContext();

    /**
     * Initialize lenses and source planes.
     *
     * @param lensRedshifts List of redshifts for each lens.
     * @param params Parameters for lenses, must be the same length as
     * lensRedshifts.
     * @param sourceRedshifts List of source plane redshifts.
     */
    int init(const std::vector<float> &lensRedshifts,
             const std::vector<std::vector<PlummerParams>> &params,
             const std::vector<float> &sourceRedshifts);

    /**
     * Set thetas for calculations.
     *
     * @param thetas List of thetas per source plane.
     */
    int setThetas(const std::vector<std::vector<Vector2D<float>>> &thetas);

    int calculatePositions(const std::vector<std::vector<float>> &masses);

    const std::vector<Vector2D<float>> &getSourcePositions(int idx) const;
};
