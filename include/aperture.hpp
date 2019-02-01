#ifndef  _APERTURE_HPP_
#define  _APERTURE_HPP_

#include "dynamic_variables.hpp"
#include "parameters.hpp"
#include "kernel/grid.hpp"

template< typename Real_t, std::size_t DGrid, std::size_t DPtc >
class Aperture {
private:
  DynamicVars<Real_t, DGrid, DPtc> _dvars;
  Params<Real_t, DGrid> _params;
  Grid<DGrid, Real_t> _grid;
  auto _comm;
  // data export

  int _timestep_begin;
  int _timestep_end;

public:
  Aperture();
  void launch();
};

#endif
