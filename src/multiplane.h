#pragma once

#include "composite.h"
#include "sources.h"
#include "util/cosmology.h"
#include <vector>

/**
 * Data for a lens plane.
 */
class PlaneData {
  public:
    CompositeLens lens;
    float redshift;

    /**
     * Construct PlaneData.
     *
     * @param l CompositeLens.
     * @param z Redshift.
     */
    PlaneData(CompositeLens &l, float z) : lens(l) { redshift = z; }

    //  __device__ PlaneData() {}

    bool operator<(const PlaneData &cmp) const {
        return redshift < cmp.redshift;
    }
};

/**
 * Multiplane lenses and source planes.
 */
class Multiplane {
  private:
    const PlaneData *__restrict__ m_plane_ptr;
    const SourcePlane *__restrict__ m_src_plane_ptr;
    PlaneData *m_plane_data;
    SourcePlane *m_src_data;
    const int m_plane_length;
    const int m_src_length;
    bool m_cuda;

  public:
    /**
     * Create new Multiplane object.
     *
     * @param plane_length Number of lens planes.
     * @param src_length Number of source planes.
     * @param plane_ptr Pointer to lens plane memory.
     * @param src_plane_ptr Pointer to source plane memory.
     */
    Multiplane(int plane_length, int src_length, PlaneData *plane_ptr,
               SourcePlane *src_plane_ptr, bool cuda = false)
        : m_plane_ptr(plane_ptr), m_src_plane_ptr(src_plane_ptr),
          m_plane_length(plane_length), m_src_length(src_length) {
        m_plane_data = plane_ptr;
        m_src_data = src_plane_ptr;
        m_cuda = cuda;
    }

    /**
     * Cleanup memory allocations.
     */
    int destroy();

    /**
     * Trace theta to source plane
     */
    __host__ __device__ uint8_t traceTheta(Vector2D<float> theta) const;
    /**
     * Trace thetas and save positions to multiple source planes
     */
    __host__ __device__ void traceTheta(Vector2D<float> theta, float *beta_x,
                                        float *beta_y,
                                        const size_t offset) const;
#ifdef __CUDACC__
    /**
     * Trace thetas and save positions to multiple source planes
     */
    __host__ __device__ void traceTheta(float2 theta, float2 *beta,
                                        const size_t offset) const;

#endif
    void traceMultiTheta(const Vector2D<float> *thetas, Vector2D<float> *betas,
                         const int length, const int plane);

    /**
     * Update lens masses (GPU only)
     *
     * @param dim LensPlane index.
     * @param i Sublens index.
     * @param masses Masses for LensPlane.
     */
    __device__ void updateLensMasses(const int dim, const int i,
                                     const double *masses) {
        m_plane_data[dim].lens.setMass(i, masses[i]);
    }

    int srcLength() const { return m_src_length; }
};

/**
 * Builder for Multiplane.
 */
class MultiplaneBuilder {
  private:
    std::vector<CompositeLensBuilder> m_builders;
	std::vector<float> m_source_z;
	std::vector<std::vector<float>> m_distances;

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
     * Add a source plane.
     *
     * @param plane SourcePlane.
     */
    void addSourcePlane(SourcePlane &plane);

    /**
     * If we aren't using source planes with points we only need the
     * redshifts for beta vectors
     */
    void setRedshifts(std::vector<float> &redshifts);

    /**
     * Get Multiplane (CPU).
     */
    Multiplane getMultiPlane();
    /**
     * Get Multiplane (GPU/CUDA).
     */
    Multiplane getCuMultiPlane();
};
