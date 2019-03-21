#ifndef  _APERTURE_HPP_
#define  _APERTURE_HPP_

#include "dynamic_variables.hpp"
#include "parameters.hpp"
#include "kernel/grid.hpp"
#include "parallel/mpi++.hpp"

template< typename Real, std::size_t DGrid, std::size_t DPtc, typename state_t >
class Aperture {
private:
  // DynamicVars<Real, DGrid, DPtc, state_t> _dvars; // TODO Field doesn't have default constructor
  Params<Real, DGrid> _params;
  // knl::Grid< Real, DGrid > _grid; // TODO Grid doesn't have default constructor
  // std::optional<mpi::Comm> _cartesian;
  // std::optional<mpi::Comm> _ensemble;
  // data export

  int _timestep_begin;
  int _timestep_end;

  struct AAA {
    void foo();
  };

public:
  Aperture();
  void launch();
};

#endif
