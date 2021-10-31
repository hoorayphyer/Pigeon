#pragma once

#include "apt/action_base.hpp"
#include "apt/grid.hpp"
#include "particle/array.hpp"
#include "particle/map.hpp"
#include "particle/properties.hpp"
#include "random/rng.hpp"

namespace field {
template <typename, int, int>
struct Field;
}

namespace dye {
template <int>
struct Ensemble;
}

namespace particle {
template <int DGrid, typename R, template <typename> class S, typename RJ>
struct Action : public apt::ActionBase<DGrid> {
  virtual ~Action(){};

  virtual Action* Clone()
      const = 0;  // covariant return types, see Modern C++ Design

  virtual void operator()(map<array<R, S>>& particles,
                          field::Field<RJ, 3, DGrid>& J,
                          std::vector<Particle<R, S>>* new_ptc_buf,
                          const map<Properties>& properties,
                          const field::Field<R, 3, DGrid>& E,
                          const field::Field<R, 3, DGrid>& B,
                          const apt::Grid<R, DGrid>& grid,
                          const dye::Ensemble<DGrid>* ens, R dt, int timestep,
                          util::Rng<R>& rng) = 0;
};
}  // namespace particle
