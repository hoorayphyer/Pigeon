#ifndef _CHECKPOINT_HPP_
#define _CHECKPOINT_HPP_

#include "dye/ensemble.hpp"
#include <string>
#include <optional>

namespace field {
  template < typename , int, int > struct Field;
}

namespace particle {
  template < typename > struct map;
  template < typename, template < typename > class > struct array;
  struct Properties;
}

// TODO other things to save, like initial conditions, shape functions? No, use gen.hpp and build a target called resume
namespace ckpt {
  template < int DGrid,
             typename Real,
             template < typename > class PtcSpecs >
  std::string save_checkpoint( std::string prefix, const int num_parts,
                               const std::optional<dye::Ensemble<DGrid>>& ens_opt,
                               int timestep,
                               const field::Field<Real, 3, DGrid>& E,
                               const field::Field<Real, 3, DGrid>& B,
                               const particle::map<particle::array<Real,PtcSpecs>>& particles,
                               const particle::map<particle::Properties>& properties,
                               const particle::map<double>& N_scat
                               );

  template < int DGrid,
             typename Real,
             template < typename > class PtcSpecs >
  int load_checkpoint( std::string dir,
                       std::optional<dye::Ensemble<DGrid>>& ens_opt,
                       const std::optional<mpi::CartComm>& cart_opt,
                       field::Field<Real, 3, DGrid>& E,
                       field::Field<Real, 3, DGrid>& B,
                       particle::map<particle::array<Real,PtcSpecs>>& particles,
                       const particle::map<particle::Properties>& properties,
                       particle::map<double>& N_scat,
                       int target_load = 0
                       );
}

#endif
