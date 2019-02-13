#ifndef  _APERTURE_HPP_
#define  _APERTURE_HPP_

#include "dynamic_variables.hpp"
#include "parameters.hpp"
#include "kernel/grid.hpp"
#include "parallel/mpi++.hpp"

template< typename Real, std::size_t DGrid, std::size_t DPtc, typename state_t >
class Aperture {
private:
  DynamicVars<Real, DGrid, DPtc, state_t> _dvars;
  Params<Real, DGrid> _params;
  knl::Grid<DGrid, Real> _grid;
  mpi::Comm _comm;
  // data export

  int _timestep_begin;
  int _timestep_end;

public:
  Aperture();
  void launch();
};

#endif
