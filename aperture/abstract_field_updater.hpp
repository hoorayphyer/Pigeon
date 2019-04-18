#ifndef  _APERTURE_ABSTRACT_FIELD_UPDATER_HPP_
#define  _APERTURE_ABSTRACT_FIELD_UPDATER_HPP_

#include "field/field.hpp"

namespace aperture {
  template < typename Real, int DGrid, typename RealJ >
  struct AbstractFieldUpdater {
    virtual void operator() ( field::Field<Real,3,DGrid>& E,
                              field::Field<Real,3,DGrid>& B,
                              const field::Field<RealJ,3,DGrid>& J,
                              Real dt, int timestep ) = 0;
  };
}

#endif
