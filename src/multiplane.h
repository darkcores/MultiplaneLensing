#pragma once

#include "composite.h"
#include "util/cosmology.h"

/**
 * Multiplane lenses and source planes.
 */
class Multiplane {
  private:
    CompositeLens *__restrict__ m_lenses;
    float *__restrict__ m_sources;
    float *__restrict__ m_dist_lenses;
    float *__restrict__ m_dist_sources;
    std::vector<int> m_dist_offsets;
    const int m_lenses_size, m_sources_size;
    const bool m_cuda;

  public:
    /**
     * Create new Multiplane object.
     *
     */
    Multiplane(CompositeLens *lenses, int lenses_size, float *sources,
               int sources_size, float *dist_lenses, float *dist_sources,
               std::vector<int> dist_offsets, bool cuda = false)
        : m_lenses_size(lenses_size), m_sources_size(sources_size),
          m_cuda(cuda) {
        m_lenses = lenses;
        m_sources = sources;
        m_dist_lenses = dist_lenses;
        m_dist_sources = dist_sources;
        m_dist_offsets = dist_offsets;
    }

    /**
     * Cleanup memory allocations.
     */
    int destroy();

/**
 * Trace theta vectors
 */
#ifdef __CUDACC__
    int traceThetas(const float2 *thetas, float2 *betas, const int n,
                    const int plane);
#else
    int traceThetas(const Vector2D<float> *thetas, Vector2D<float> *betas,
                    const int n, const int plane);
#endif
};

/**
 * Builder for Multiplane.
 */
class MultiplaneBuilder {
  private:
    std::vector<CompositeLensBuilder> m_builders;
    std::vector<float> m_source_z;
    std::vector<float> m_dists_lenses;
    std::vector<float> m_dists_sources;
    std::vector<int> m_dist_offsets;

    const Cosmology m_cosm;

    void prepare();

  public:
    /**
     * New Multiplane builder.
     *
     * @param cosm Cosmology settings.
     */
    MultiplaneBuilder(const Cosmology cosm) : m_cosm(cosm) {}

    /**
     * Add a lens plane.
     *
     * @param lensbuilder Builder for CompositeLens.
     */
    void addPlane(CompositeLensBuilder &lensbuilder);

    /**
     * If we aren't using source planes with points we only need the
     * redshifts for beta vectors
     */
    void setRedshifts(std::vector<float> &redshifts) { m_source_z = redshifts; }

    /**
     * Get Multiplane (CPU).
     */
    Multiplane getMultiPlane();
    /**
     * Get Multiplane (GPU/CUDA).
     */
    Multiplane getCuMultiPlane();
};
