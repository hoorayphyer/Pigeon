#ifndef  _ABSTRACT_FIELD_UPDATER_HPP_
#define  _ABSTRACT_FIELD_UPDATER_HPP_

#include "field/field.hpp"

namespace pic {
  template < typename Real, int DGrid, typename RealJ >
  struct AbstractFieldUpdater {
    virtual void operator() ( field::Field<Real,3,DGrid>& E,
                              field::Field<Real,3,DGrid>& B,
                              const field::Field<RealJ,3,DGrid>& J,
                              Real dt, int timestep ) = 0;
  };
}

#endif
