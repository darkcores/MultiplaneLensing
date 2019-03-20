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
 */
class MultiPlaneContext {
  private:
    const double m_angularUnit;
    const Cosmology m_cosmology;
    // Device thetas
    float *m_theta_x, *m_theta_y;
    size_t m_theta_len;
    // Device betas, length is m_theta_len * number of source planes
    float *m_beta_x, *m_beta_y;
    Multiplane *m_multiplane;

    CompositeLensBuilder buildLens(const double Dd,
                                   const std::vector<PlummerParams> &params);

  public:
    MultiPlaneContext(const double angularUnit, const Cosmology cosmology);
    ~MultiPlaneContext();

    int init(const std::vector<float> &lensRedshifts,
             const std::vector<std::vector<PlummerParams>> &params,
             const std::vector<float> &sourceRedshifts);

    int setThetas(const std::vector<Vector2D<float>> &thetas);

    int calculatePositions(const std::vector<std::vector<double>> &masses);

    const std::vector<Vector2D<float>> getSourcePositions(int idx) const;
};
