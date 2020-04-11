#ifndef  _PARTICLE_ANNIHILATION_HPP_
#define  _PARTICLE_ANNIHILATION_HPP_

#include "particle/array.hpp"
#include "field/field.hpp"
#include "apt/grid.hpp"

namespace mpi { struct Comm; }

namespace particle {
  template < int DGrid, typename R, template < typename > class S, typename ShapeF, typename RJ >
  void annihilate( array<R,S>& el, array<R,S>& po,
                   field::Field<RJ,3,DGrid>& J,
                   R charge_el, R charge_po,
                   const apt::Grid< R, DGrid >& grid,
                   const mpi::Comm& intra,
                   R dt, const ShapeF&,
                   R(*policy)(R num_electron_in_a_cell, R num_positron_in_a_cell),
                   const apt::array<apt::array<R,2>,DGrid>& bounds );
}

#endif
