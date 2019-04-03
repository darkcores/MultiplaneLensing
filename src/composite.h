#pragma once

#include "plummer.h"
#include <vector>

class CompositeLens {
  private:
    Plummer *__restrict__ m_lenses;
    const int m_lenses_size;
    const bool m_cuda;

  public:
    CompositeLens(Plummer *lenses, const int lenses_size,
                  const bool cuda = true)
        : m_lenses(lenses), m_lenses_size(lenses_size), m_cuda(cuda) {}

    int destroy();
};

/**
 * Builder for CompositeLens. Please beware currently also manager
 * memory for CompositeLens, will change in the future. For now, do
 * not delete until CompositeLens is no longer needed.
 */
class CompositeLensBuilder {
  private:
    std::vector<Plummer> m_lenses;

    float m_redshift;
    float m_scale;

  public:
    /**
     * New CompositeLensBuilder.
     */
    CompositeLensBuilder(const float scale, const float redshift) {
        m_scale = scale;
        m_redshift = redshift;
    }

    float redshift() const { return m_redshift; }

    /**
     * Add plummer sublens.
     *
     * @param lens Plummer sublens.
     * @param position Plummer sublens position.
     */
    void addLens(const Plummer &lens);

    /**
     * Get cuda lens. Has internal pointer to cuda memory.
     */
    CompositeLens getCuLens();
    /**
     * Get lens. Has internal pointer to host memory.
     */
    CompositeLens getLens();
};
