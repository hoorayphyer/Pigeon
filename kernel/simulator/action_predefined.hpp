#pragma once

#include "apt/grid.hpp"
#include "simulator/action.hpp"

namespace particle {

template <int DGrid, typename R, template <typename> class S, typename RJ>
class Migrator
    : public pic::ActionWithSetters<Migrator<DGrid, R, S, RJ>, DGrid,
                                    pic::ParticleAction<DGrid, R, S, RJ>> {
 public:
  using Bundle_t = pic::ParticleAction<DGrid, R, S, RJ>::Bundle_t;

  auto& set_supergrid(const apt::Grid<R, DGrid>& supergrid) noexcept {
    m_supergrid = supergrid;
    return *this;
  }

  void operator()(const Bundle_t& bundle) const override {
    impl(bundle.particles, bundle.new_ptc_buf, bundle.properties, bundle.grid,
         bundle.ens, bundle.timestep);
  }

 private:
  apt::Grid<R, DGrid> m_supergrid;

  // FIXME migrate need more memory check
  void impl(map<array<R, S>>& particles,
            std::vector<Particle<R, S>>& new_ptc_buf,
            const map<Properties>& properties, const apt::Grid<R, DGrid>& grid,
            const dye::Ensemble<DGrid>& ens, int timestep) const;
};

}  // namespace particle
