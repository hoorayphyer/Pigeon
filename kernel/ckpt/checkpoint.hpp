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
}

// TODO other things to save, like initial conditions, shape functions? No, use gen.hpp and build a target called resume
namespace ckpt {
  template < int DGrid,
             typename Real,
             template < typename > class PtcSpecs >
  void save_checkpoint( std::string prefix, const int num_parts,
                        const std::optional<dye::Ensemble<DGrid>>& ens_opt,
                        int timestep,
                        const field::Field<Real, 3, DGrid>& E,
                        const field::Field<Real, 3, DGrid>& B,
                        const particle::map<particle::array<Real,PtcSpecs>>& particles
                        );

  // void save_tracing();
}

#endif
