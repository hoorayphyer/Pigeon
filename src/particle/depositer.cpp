#include "particle/depositer.hpp"
#include "apt/numeric.hpp"
#include "apt/block.hpp"
#include "apt/vec.hpp"

namespace esirkepov {
  template < typename T >
  inline T Wesir_2D( T sx0, T sx1, T sy0, T sy1 ) noexcept {
    return ( ( 2 * sx1 + sx0 ) * sy1 + ( sx1 + 2 * sx0 ) * sy0 ) / 6.0;
  }

  template < typename T >
  inline T Wesir( T sx0, T sx1, T sy0, T sy1, T sz0, T sz1 ) noexcept {
    return (sx1 - sx0) * Wesir_2D(sy0, sy1, sz0, sz1);
  }

  // TODO make sure ShapeF can be passed in as constexpr. This may be used to optimize the ever-checking away.
  template < typename T, int DField, int DGrid, typename ShapeF >
  constexpr std::array<T,DField> calcW( const apt::Index<DGrid>& index,
                                        const std::array<T, DGrid>& sep0_b,
                                        const std::array<T, DGrid>& sep1_b,
                                        const ShapeF& shapef ) noexcept {
    static_assert( DField == 3 );
    T W[DField]{};

    // TODO use normalization?? er ci xing zheng ze hua
    T s0[DGrid]{};
    T s1[DGrid]{};

    auto f_calc_s0s1
      = [&shapef] ( auto& s0_i, auto& s1_i, auto i, auto sep0_i, auto sep1_i ) {
          s0_i = shapef( i + sep0_i );
          s1_i = shapef( i + sep1_i );
        };

    apt::foreach<0,DGrid>( f_calc_s0s1, s0, s1, index, sep0_b, sep1_b );

    if constexpr ( DGrid == 2 ) {
        W[0] = Wesir( s0[0], s1[0], s0[1], s1[1], 1.0, 1.0 );
        W[1] = Wesir( s0[1], s1[1], 1.0, 1.0, s0[0], s1[0] );
        W[2] = Wesir( 0.0, 1.0, s0[0], s1[0], s0[1], s1[1] );
      } else {
      W[0] = Wesir( s0[0], s1[0], s0[1], s1[1], s0[2], s1[2] );
      W[1] = Wesir( s0[1], s1[1], s0[2], s1[2], s0[0], s1[0] );
      W[2] = Wesir( s0[2], s1[2], s0[0], s1[0], s0[1], s1[1] );
    }

    return {{ W[0], W[1], W[2] }};
  }

}

namespace particle {

  // RATIONALE: assume the offset = 0.5 in the dimension. The native grid is one offset by 0.5 with respect to the original one.
  // - q1 is relative position in the native grid. q1 starts at 0.0
  // - the contributing cells in the native grid are [ int(q1 - sf.r) + 1,  int(q1 + sf.r) + 1 )
  // - index_original = index_native + ( q1 - int(q1) >= 0.5 )
  // Now we have q0 and q1, the final range of contributing cells is the union of individual ones
  template < typename Vec_q1_abs, typename Vec_dq_abs, typename Grid, typename ShapeF >
  constexpr auto prep( const Vec_q1_abs& q1_abs, const Vec_dq_abs& dq_abs, const Grid& grid, const ShapeF& ) {
    constexpr int DGrid = apt::ndim_v<Grid>;
    using T = apt::most_precise_t<apt::element_t<Vec_q1_abs>, apt::element_t<Vec_dq_abs>>;
    apt::Index<DGrid> I_b;
    apt::Index<DGrid> extent;
    std::array<T, DGrid> sep0_b; // separation 0 is defined to be i_cell_beginning - q0_ptc
    std::array<T, DGrid> sep1_b;

    apt::foreach<0, DGrid>
      ( [] ( auto& i_b, auto& ext, auto& sp0, auto& sp1,
             auto qmin, auto qmax, const auto& g1d )
        noexcept {
          // initially, qmin = q1_abs, qmax = dq_abs.
          qmax /= g1d.delta();
          qmin = ( qmin - g1d.lower() ) / g1d.delta() - 0.5;
          // now, qmin = q1_rel, qmax = dq_rel.
          sp1 = - qmin;
          sp0 = sp1 + qmax;
          qmin -= (( qmax > 0 ? qmax : 0.0 ) + ShapeF::support / 2.0);
          qmax = qmin + qmax * ( (qmax > 0.0) - ( qmax < 0.0) ) + ShapeF::support;
          i_b = int( qmin ) + 1 ;
          ext = int( qmax ) + 1 - i_b;
          sp1 += i_b;
          sp0 += i_b;
        },
        I_b, extent, sep0_b, sep1_b,
        q1_abs, dq_abs, grid );

    return std::make_tuple( I_b, extent, sep0_b, sep1_b );
  }

}

namespace particle {

  template < typename Field,
             typename Ptc,
             typename Vec,
             typename ShapeF >
  void depositWJ ( Field& WJ,
                   const PtcExpression<Ptc>& ptc,
                   const apt::VecExpression<Vec>& dq_ptc,
                   const ShapeF& shapef ) {
    constexpr int DGrid = apt::ndim_v<decltype(WJ.mesh().bulk())>;
    static_assert( DGrid > 1 && DGrid < 4 );


    const auto[I_b, extent, sep0_b, sep1_b] = prep( ptc.q(), dq_ptc, WJ.mesh().bulk(), shapef );

    for ( const auto& I : apt::Block(extent) ) {
      auto W = esirkepov::calcW( I, sep0_b, sep1_b, shapef );
      apt::foreach<0, DGrid>
        ( [&]( auto& wj, auto w ) { wj(I_b + I) += w; } , WJ, W );

      // FIXME fix the following in DGrid == 2
      // Calling deposition after pusher.calculateDisplacement implies
      // that p_tmp, which is at n+0.5 time step, is based with respect
      // to x^(n+1). However, the x used here is actually x^(n). One way
      // to improve this is obviously calling updatePos before deposition
      // and accordingly change expressions for calculating shapefunctions.
      // FIXME: But, where is J based? Does one really need rebasing momentum?
      if constexpr ( DGrid == 2 )
                     std::get<2>(WJ)( I_b + I ) += std::get<2>(W) * std::get<2>(ptc.p()) / std::sqrt( 1.0 + apt::sqabs(ptc.p()) );

    }

    return;
  }

}
