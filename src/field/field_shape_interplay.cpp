#include "field/field_shape_interplay.hpp"
#include "apt/block.hpp"
#include "apt/vec.hpp"
#include "apt/type_traits.hpp"
#include <tuple>

namespace esirkepov {
  template < typename T >
  inline T Wesir_2D( T sx0, T sx1, T sy0, T sy1 ) noexcept {
    return ( ( 2 * sx1 + sx0 ) * sy1 + ( sx1 + 2 * sx0 ) * sy0 ) / 6.0;
  }

  template < typename T >
  inline T Wesir( T sx0, T sx1, T sy0, T sy1, T sz0, T sz1 ) noexcept {
    return (sx1 - sx0) * Wesir_2D(sy0, sy1, sz0, sz1);
  }

  template < typename T, int DGrid, typename ShapeF >
  constexpr std::array<T,3> calcW( const apt::Index<DGrid>& index,
                                   const std::array<T, DGrid>& sep0_b,
                                   const std::array<T, DGrid>& sep1_b,
                                   const ShapeF& shapef ) noexcept {
    T W[3]{};

    // TODOL use normalization?? er ci xing zheng ze hua
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

namespace field {

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

namespace field {

  template < typename Field, typename Vec_q, typename Vec_dq, typename ShapeF >
  void depositWJ ( Field& WJ,
                   const apt::VecExpression<Vec_q>& q1_abs,
                   const apt::VecExpression<Vec_dq>& dq_abs,
                   const ShapeF& shapef ) {
    constexpr int DGrid = apt::ndim_v<decltype(WJ.mesh().bulk())>;
    constexpr int DField = apt::ndim_v<Field>;
    static_assert( DField == 3 );
    static_assert( DGrid > 1 && DGrid < 4 );

    const auto[I_b, extent, sep0_b, sep1_b] = prep( q1_abs, dq_abs, WJ.mesh().bulk(), shapef );

    for ( const auto& I : apt::Block(extent) ) {
      auto W = esirkepov::calcW( I, sep0_b, sep1_b, shapef );
      if constexpr ( DGrid == 2 ) std::get<2>(W) *= std::get<2>(dq_abs); // see NOTE below
      apt::foreach<0, DField> // NOTE it is DField here, not DGrid
        ( [&]( auto& wj, auto w ) { wj(I_b + I) += w; } , WJ, W );
    }

    // NOTE in eq.36 in Esirkepov, V_z is really del_z / dt, where del_z should
    // be interpretted as the COORDINATE displacement of the particle and dt is
    // the timestep. It is wrong to use the physical velcocity v_ptc_z here,
    // which is what we did before in LogSpherical where there was always a
    // confusion about where to base the v_ptc_z.

    return;
  }

}

namespace field {
  template < typename T, int DField, int DGrid, typename Vec_q, typename ShapeF >
  apt::Vec<T, DField> interpolate ( const Field<T,DField,DGrid>& field,
                                    const apt::VecExpression<Vec_q>& q_abs,
                                    const ShapeF& shapef ) {
    apt::Vec<T, DField> result;
    const auto& grid = field.mesh().bulk();
    constexpr auto supp = ShapeF::support;
    std::array<T, DGrid> loc {}; // location at which to interpolate field

    apt::foreach<0,DGrid>
      ( []( auto& l, auto q, const auto& g ){
          l = ( q - g.lower() ) / g.delta();
        }, loc, q_abs, grid );

    apt::Index<DGrid> I_b {};

    apt::foreach<0,DField>
      ( [&] ( auto& res, const auto& comp ) {
          res = 0.0;
          std::array<T,DGrid> sep_b {};

          // correct loc by offset
          apt::foreach<0,DGrid>
            ( []( auto& l, auto ofs ){
                l -= ofs;
              }, loc, comp.offset );

          apt::foreach<0,DGrid>
            ( []( auto& i_b, auto& s_b, auto l ){
                i_b = int(l - supp / 2.0) + 1;
                s_b = i_b - l;
              }, I_b, sep_b, loc );

          T wgt [DGrid] {};

          auto calc_wgt
            = [&]( auto& wgt, const auto& I ) {
                wgt[0] = shapef( I[0] + sep_b[0] );
                if constexpr ( DGrid > 1 ) {
                    if( I[0] != 0 ) return;
                    wgt[1] = shapef( I[1] + sep_b[1] );
                  }
                if constexpr ( DGrid > 2 ) {
                    if( I[1] != 0 ) return;
                    wgt[2] = shapef( I[2] + sep_b[2] );
                  }
              };

          for ( const auto& I : apt::Block( { supp, supp, supp } ) ) {
            calc_wgt( wgt, I );
            res += comp( I_b + I ) * wgt[0] * wgt[1] * wgt[2];
          }

          // restore loc
          apt::foreach<0,DGrid>
            ( []( auto& l, auto ofs ){
                l += ofs;
              }, loc, comp.offset );

        }, result, field );


    return result;
  }

}
