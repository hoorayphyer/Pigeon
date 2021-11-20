#pragma once

#include <unordered_set>

#include "apt/grid.hpp"
#include "particle/array.hpp"
#include "particle/map.hpp"
#include "particle/properties.hpp"
#include "particle/species_predef.hpp"
#include "random/rng.hpp"

namespace particle {
template <int DGrid, typename R, template <typename> class S, typename ShapeF,
          typename RJ>
class Updater {
 public:
  using UpdateQ_t = apt::array<R, S<R>::Dim> (*)(
      typename array<R, S>::particle_type::vec_type& x,
      typename array<R, S>::particle_type::vec_type& p, R dt, bool is_massive);

  Updater& set_update_q(UpdateQ_t update_q) {
    _update_q = update_q;
    return *this;
  }

  Updater& set_ignore_current(std::unordered_set<species> sps) {
    _ignore_current = std::move(sps);
    return *this;
  }

  void operator()(map<array<R, S>>& particles, field::Field<RJ, 3, DGrid>& J,
                  std::vector<Particle<R, S>>* new_ptc_buf,
                  const map<Properties>& properties,
                  const field::Field<R, 3, DGrid>& E,
                  const field::Field<R, 3, DGrid>& B,
                  const apt::Grid<R, DGrid>& grid, R dt, int timestep,
                  util::Rng<R>& rng) const;

 private:
  UpdateQ_t _update_q;
  std::unordered_set<species> _ignore_current{};
};
}  // namespace particle
