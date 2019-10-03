#ifndef _PIC_HPP_
#define _PIC_HPP_

#include "particle/shapef.hpp"
#include "apt/type_traits.hpp"

namespace particle {
  // TODOL users may forget to sync value and state. Add another layer then
  template < typename T >
  struct Specs {
    using value_type = T;
    static constexpr int Dim = 3;
    using state_type = apt::copy_cvref_t<T, long long>;

    static_assert( 8 * sizeof( state_type ) >= 64 );
  };
}

namespace pic {
  constexpr int DGrid = 2;
  using real_t = float;
  using real_j_t = float;
  using real_export_t = float;

  using ShapeF = particle::shapef_t<particle::shape::Cloud_In_Cell>;
};


#endif
