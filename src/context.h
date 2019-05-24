#pragma once

#include <vector>
// Cosmology calculations are needed anyways
#include "util/cosmology.h"
// Same for points/Vector2D
#include "util/vector2d.h"

// Older GCC versions don't work otherwise
#ifndef size_t
typedef std::size_t size_t;
#endif

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
 * Return a how many cuda devices are available.
 * @param print Prints device info for each device to stdout
 */
int cudaDevices(bool print = false);

/**
 * Context for multiplane plummer calculations.
 *
 * For sample use cases see api_tests.cu in tests/ or example.cpp in
 * example/ for a complete example program.
 */
class MultiPlaneContext {
  private:
	const int m_device;
    const double m_angularUnit;
    const Cosmology m_cosmology;
    Vector2D<float> *m_theta;
    size_t m_theta_len, m_source_len;
    std::vector<size_t> m_theta_count, m_lens_count;
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
	 * @param device CUDA device id, used for multi GPU. This class is 
	 * synchronous, so using a thread per device is recommended.
     */
    MultiPlaneContext(const double angularUnit, const Cosmology cosmology,
                      const int device = 0);
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
    // int setThetas(const std::vector<std::vector<Vector2D<float>>> &thetas);
    int setThetas(const std::vector<std::vector<Vector2D<float>>> &thetas);

    /**
     * Calculation beta positions with lens masses.
     */
    int calculatePositions(
        const std::vector<std::vector<float>> &masses,
        const std::vector<float> &mass_sheet = std::vector<float>(0));

    int
    calculatePositionsBenchmark(const std::vector<std::vector<float>> &masses,
                                float &millis, int nruns = 7);

    const std::vector<Vector2D<float>> &getSourcePositions(int idx) const;
};
