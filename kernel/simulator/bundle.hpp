#pragma once
#include <vector>

#include "apt/grid.hpp"
#include "dye/ensemble.hpp"
#include "field/field.hpp"
#include "mpipp/mpi++.hpp"
#include "particle/array.hpp"
#include "particle/map.hpp"
#include "particle/particle.hpp"
#include "particle/properties.hpp"
#include "random/rng.hpp"

namespace pic {

template <int DGrid, typename R, typename RJ>
struct FieldBundle {
  field::Field<R, 3, DGrid>& E;
  field::Field<R, 3, DGrid>& B;
  field::Field<RJ, 3, DGrid>& J;
  const apt::Grid<R, DGrid>& grid;
  const mpi::CartComm& cart;
  int timestep;
  R dt;
};

template <int DGrid, typename R, template <typename> class S, typename RJ>
struct ParticleBundle {
  particle::map<particle::array<R, S>>& particles;
  field::Field<RJ, 3, DGrid>& J;
  std::vector<particle::Particle<R, S>>& new_ptc_buf;
  const particle::map<particle::Properties>& properties;
  const field::Field<R, 3, DGrid>& E;
  const field::Field<R, 3, DGrid>& B;
  const apt::Grid<R, DGrid>& grid;
  const dye::Ensemble<DGrid>& ens;
  R dt;
  int timestep;
  util::Rng<R>& rng;
};

template <int DGrid, typename R, template <typename> class S, typename RJ>
struct ExportBundle {
  const particle::map<particle::array<R, S>>& particles;
  const particle::map<particle::Properties>& properties;
  const field::Field<R, 3, DGrid>& E;
  const field::Field<R, 3, DGrid>& B;
  const field::Field<RJ, 3, DGrid>& J;
  const apt::Grid<R, DGrid>& grid;
  const std::optional<mpi::CartComm>& cart_opt;
  const dye::Ensemble<DGrid>& ens;
  R dt;
  int timestep;
};

// RATIONALE
//
// InitialConditionBundle and PostResumeBundle are mutually exclusive. That
// is, depending on the presence of a resume dir, only one of them will need
// to be executed. From the program's point of view, we only need to maintain
// one Bundle class. But this way user will have to modify their code just to
// make sure to use only the correct actions, which is too error prone.
// Therefore, we keep these separate

template <int DGrid, typename R, template <typename> class S, typename RJ>
struct InitialConditionBundle {
  field::Field<R, 3, DGrid>& E;
  field::Field<R, 3, DGrid>& B;
  field::Field<RJ, 3, DGrid>& J;
  particle::map<particle::array<R, S>>& particles;
  const particle::map<particle::Properties>& properties;
  const apt::Grid<R, DGrid>& grid;
};

template <int DGrid, typename R, template <typename> class S, typename RJ>
struct PostResumeBundle {
  field::Field<R, 3, DGrid>& E;
  field::Field<R, 3, DGrid>& B;
  field::Field<RJ, 3, DGrid>& J;
  particle::map<particle::array<R, S>>& particles;
  const particle::map<particle::Properties>& properties;
  const apt::Grid<R, DGrid>& grid;
  const std::optional<dye::Ensemble<DGrid>>& ens_opt;
  int resumed_timestep;
  std::string this_run_dir;
};

}  // namespace pic
