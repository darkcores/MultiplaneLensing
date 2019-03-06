#pragma once

#include "composite.h"
#include "sources.h"
#include "util/cosmology.h"
#include <vector>

#include <cuda_runtime_api.h>

class PlaneData {
  public:
    CompositeLens lens;
    float redshift;

    PlaneData(CompositeLens &l, float z) : lens(l) { redshift = z; }

    __device__ PlaneData() {}

    bool operator<(const PlaneData &cmp) const {
        return redshift < cmp.redshift;
    }
};

class Multiplane {
  private:
    const PlaneData *__restrict__ m_plane_ptr;
    const SourcePlane *__restrict__ m_src_plane_ptr;
    const int m_plane_length;
    const int m_src_length;

  public:
    Multiplane()
        : m_plane_ptr(nullptr), m_src_plane_ptr(nullptr), m_plane_length(0),
          m_src_length(0) {}
    Multiplane(int plane_length, int src_length, PlaneData *plane_ptr,
               SourcePlane *src_plane_ptr)
        : m_plane_ptr(plane_ptr), m_src_plane_ptr(src_plane_ptr),
          m_plane_length(plane_length), m_src_length(src_length) {}

    void getImage() const;

    /**
     * Trace theta to source plane
     */
    uint8_t traceTheta(Vector2D<float> theta) const;
};

class MultiplaneBuilder {
  private:
	std::vector<CompositeLensBuilder *> m_builders;
	std::vector<PlaneData> m_data;
	std::vector<SourcePlane> m_src_data;

    PlaneData *plane_ptr;
    SourcePlane *src_ptr;

	const Cosmology m_cosm;

    void prepare();

  public:
    MultiplaneBuilder(const Cosmology cosm) : m_cosm(cosm) {}

    void addPlane(CompositeLensBuilder *lensbuilder);
    void addSourcePlane(SourcePlane &plane);

    Multiplane getMultiPlane();
    Multiplane getCuMultiPlane();
	void cuFree();
};
