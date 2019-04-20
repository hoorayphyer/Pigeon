#ifndef _PIC_HPP_
#define _PIC_HPP_

#include "kernel/shapef.hpp"
#include "kernel/curvilinear.hpp"
#include "apt/type_traits.hpp"

namespace particle {
  // TODOL users may forget to sync value and state. Add another layer then
  template < typename T >
  struct Specs {
    using value_type = T;
    static constexpr int Dim = 3;
    using state_type = apt::copy_cvref_t<T,unsigned long long>;

    static_assert( 8 * sizeof( state_type ) >= 64 );
  };
}

namespace pic {
  using real_t = double;

  constexpr int DGrid = 2;

  using ShapeF = knl::shapef_t<knl::shape::Cloud_In_Cell>;

  using Metric = knl::coord<knl::coordsys::Cartesian>;

  using real_j_t = long double;

};


#endif
