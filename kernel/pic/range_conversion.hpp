#ifndef _PIC_RANGE_CONVERSION_HPP_
#define _PIC_RANGE_CONVERSION_HPP_

namespace pic {
  // TODO check gtl boundary on extent.
  // NOTE range is assumed to be [,), and range[1] >= range[0]
  template < typename Real >
  apt::pair< int > gtl ( const apt::pair<Real>& range,
                          const mani::Grid1D<Real>& localgrid ) noexcept {
    auto to_grid_index =
      [&localgrid] ( auto absc ) noexcept {
        return ( absc - localgrid.lower() ) / localgrid.delta();
      };
    int lb = to_grid_index( range[LFT] );
    int ub = to_grid_index( range[RGT] );

    if ( ub <= 0 || lb >= localgrid.dim() ) {
      // range not applicable on current local patch
      return {0,0};
    } else {
      lb = std::max<int>( 0, lb );
      ub = std::min<int>( ub, localgrid.dim() );
      return { lb, ub - lb };
    }
  }

  template < typename Real >
  apt::pair< int > gtl ( int Ib_global, int extent_global,
                         const mani::Grid1D<Real>& supergrid,
                         const mani::Grid1D<Real>& localgrid ) noexcept {
    // lb is the Ib_global with respect to the localgrid
    int lb = Ib_global - static_cast<int>( ( localgrid.lower() - supergrid.lower() ) / localgrid.delta() + 0.5 );
    int ub = lb + extent_global;

    if ( ub <= 0 || lb >= localgrid.dim() ) {
      // range not applicable on current local patch
      return {0,0};
    } else {
      lb = std::max<int>( 0, lb );
      ub = std::min<int>( ub, localgrid.dim() );
      return { lb, ub - lb };
    }
  }

}

#endif
