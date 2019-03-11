#pragma once

#include "composite.h"
#include "sources.h"
#include "util/cosmology.h"
#include <vector>

#include <cuda_runtime_api.h>

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
    const int m_plane_length;
    const int m_src_length;

  public:
	/**
	 * New empty source plane, may be removed later.
	 */
	[[deprecated]]
    Multiplane()
        : m_plane_ptr(nullptr), m_src_plane_ptr(nullptr), m_plane_length(0),
          m_src_length(0) {}
	/**
	 * Create new Multiplane object.
	 *
	 * @param plane_length Number of lens planes.
	 * @param src_length Number of source planes.
	 * @param plane_ptr Pointer to lens plane memory.
	 * @param src_plane_ptr Pointer to source plane memory.
	 */
    Multiplane(int plane_length, int src_length, PlaneData *plane_ptr,
               SourcePlane *src_plane_ptr)
        : m_plane_ptr(plane_ptr), m_src_plane_ptr(src_plane_ptr),
          m_plane_length(plane_length), m_src_length(src_length) {
		m_plane_data = plane_ptr;
	}

    /**
     * Trace theta to source plane
     */
    __host__ __device__ uint8_t traceTheta(Vector2D<float> theta) const;

	/**
	 * Update lens masses (GPU only)
	 *
	 * @param dim LensPlane index.
	 * @param i Sublens index.
	 * @param masses Masses for LensPlane.
	 */
	__device__ void updateLensMasses(const int dim, const int i, const float *masses) {
		m_plane_data[dim].lens.setMass(i, masses[i]);
	}
};

/**
 * Builder for Multiplane.
 */
class MultiplaneBuilder {
  private:
	std::vector<CompositeLensBuilder *> m_builders;
	std::vector<PlaneData> m_data;
	std::vector<SourcePlane> m_src_data;

	bool cuda = false;
    PlaneData *plane_ptr;
    SourcePlane *src_ptr;

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
	 * Cleanup, also cleans up GPU resources. Issue see CompositeLens.
	 */
	~MultiplaneBuilder() {
		if (cuda)
			cuFree();
	}

	/**
	 * Add a lens plane.
	 * 
	 * @param lensbuilder Builder for CompositeLens.
	 */
    void addPlane(CompositeLensBuilder *lensbuilder);
	/**
	 * Add a source plane.
	 *
	 * @param plane SourcePlane.
	 */
    void addSourcePlane(SourcePlane &plane);

	/**
	 * Get Multiplane (CPU).
	 */
    Multiplane getMultiPlane();
	/**
	 * Get Multiplane (GPU/CUDA).
	 */
    Multiplane getCuMultiPlane();
	/**
	 * Free cuda allocations.
	 */
	void cuFree();
};
