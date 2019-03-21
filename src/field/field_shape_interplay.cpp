#include "field/field_shape_interplay.hpp"
#include "field/field.hpp"
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
  constexpr apt::array<T,3> calcW( const apt::Index<DGrid>& index,
                                   const apt::array<T, DGrid>& sep0_b,
                                   const apt::array<T, DGrid>& sep1_b,
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
  // - q is relative position in the native grid. q starts at 0.0
  // - the contributing cells in the native grid are [ int(q - sf.r) + 1,  int(q + sf.r) + 1 )
  // - index_original = index_native + ( q - int(q) >= 0.5 )
  // Now we have q0 and q1, the final range of contributing cells is the union of individual ones
  template < typename Vec_q_abs, typename Vec_dq_abs, typename Grid, typename ShapeF >
  constexpr auto deposit_dJ_prep( const Vec_q_abs& q0_abs, const Vec_dq_abs& dq_abs, const Grid& grid, const ShapeF& ) {
    constexpr int DGrid = apt::ndim_v<Grid>;
    using T = apt::most_precise_t<apt::element_t<Vec_q_abs>, apt::element_t<Vec_dq_abs>>;
    apt::Index<DGrid> I_b;
    apt::Index<DGrid> extent;
    apt::array<T, DGrid> sep0_b; // separation 0 is defined to be i_cell_beginning - q0_ptc
    apt::array<T, DGrid> sep1_b;

    apt::foreach<0, DGrid>
      ( [] ( auto& i_b, auto& ext, auto& sp0, auto& sp1,
             auto qmin, auto qmax, const auto& g1d )
        noexcept {
          // initially, qmin = q0_abs, qmax = dq_abs.
          qmax /= g1d.delta();
          qmin = ( qmin - g1d.lower() ) / g1d.delta() - 0.5;
          // now, qmin = q0_native, qmax = dq_native.
          sp0 = - qmin;
          sp1 = sp0 - qmax;
          qmin += ( ( qmax < 0 ? qmax : 0.0 ) - ShapeF::support / 2.0 );
          qmax = qmin + qmax * ( (qmax > 0.0) - ( qmax < 0.0) ) + ShapeF::support;
          i_b = int( qmin ) + 1 ;
          ext = int( qmax ) + 1 - i_b;
          sp0 += i_b;
          sp1 += i_b;
        },
        I_b, extent, sep0_b, sep1_b,
        q0_abs, dq_abs, grid );

    return std::make_tuple( I_b, extent, sep0_b, sep1_b );
  }

}

namespace field {

  template < typename Field, typename T, typename Vec_q, typename Vec_dq, typename ShapeF >
  void deposit_dJ ( Field& dJ, T charge,
                   const apt::VecExpression<Vec_q>& q0_abs,
                   const apt::VecExpression<Vec_dq>& dq_abs,
                   const ShapeF& shapef ) {
    constexpr int DGrid = apt::ndim_v<decltype(dJ.mesh().bulk())>;
    constexpr int DField = apt::ndim_v<Field>;
    static_assert( DField == 3 );
    static_assert( DGrid > 1 && DGrid < 4 );

    const auto[I_b, extent, sep0_b, sep1_b] = deposit_dJ_prep( q0_abs, dq_abs, dJ.mesh().bulk(), shapef );

    for ( const auto& I : apt::Block(extent) ) {
      auto W = esirkepov::calcW( I, sep0_b, sep1_b, shapef );
      if constexpr ( DGrid == 2 ) W[2] *= dq_abs[2]; // see NOTE below
      apt::foreach<0, DField> // NOTE it is DField here, not DGrid
        ( [&]( auto& dj, auto w ) { dj(I_b + I) += w * charge; } , dJ, W );
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
  namespace impl {
    template < typename T, int DGrid, typename ShapeF >
    struct WeightFinder {
    private:
      apt::Index<DGrid> _I_b {};
      apt::array<T,DGrid> _sep_b {};
      apt::array<T,DGrid> _wgt {};
      const ShapeF& _shapef;
    public:
      // NOTE loc is the relative index
      constexpr WeightFinder( const apt::array<T, DGrid>& loc,
                              const apt::array< offset_t, DGrid >& offset,
                              const ShapeF& shapef ) noexcept
        : _shapef( shapef ) {

        apt::foreach<0,DGrid>
          ( []( auto& i_b, auto& s_b, auto l, const auto& ofs ){
              l -= ofs; // now l is the native grid index
              i_b = int(l - ShapeF::support / 2.0) + 1; // i_b is native
              s_b = i_b - l;
              i_b += (( ofs == MIDWAY ) && ( l - static_cast<int>(l) >= ofs )); // i_b now is with respect to original grid
            }, _I_b, _sep_b, loc, offset );
      }

      constexpr const auto& I_b() const noexcept { return _I_b; }

      constexpr T weight ( const apt::Index<DGrid>& I ) noexcept{
        T result = 0;

        auto update_wgt = // alternative to nested loops
          [&]( const auto& I ) {
            static_assert( DGrid > 0 );
            _wgt[0] = _shapef( I[0] + _sep_b[0] );
            if constexpr ( DGrid > 1 ) {
                if( I[0] != 0 ) return; // no need to recalculate wgt[1]
                _wgt[1] = _shapef( I[1] + _sep_b[1] );
              }
            if constexpr ( DGrid > 2 ) {
                if( I[1] != 0 ) return;
                _wgt[2] = _shapef( I[2] + _sep_b[2] );
              }
          };
        update_wgt();
        apt::foreach<0,DGrid>
          ( [&result]( auto w ){ result *= w; }, _wgt );
        return result;
      };
    };
  }

  template < typename T, int DField, int DGrid, typename LocType, typename ShapeF >
  apt::Vec<T, DField> interpolate ( const Field<T,DField,DGrid>& field,
                                    const LocType& q_abs,
                                    const ShapeF& shapef ) {
    apt::Vec<T, DField> result;
    const auto& grid = field.mesh().bulk();
    constexpr auto supp = ShapeF::support;
    apt::array<T, DGrid> loc {}; // location at which to interpolate field

    apt::foreach<0,DGrid>
      ( []( auto& l, auto q, const auto& g ){
          l = ( q - g.lower() ) / g.delta();
        }, loc, q_abs, grid );

    apt::foreach<0,DField>
      ( [&] ( auto& res, const auto& comp ) {
          res = 0.0;
          auto wf = WeightFinder( loc, comp.offset, shapef );

          for ( const auto& I : apt::Block( apt::Index<DGrid>{ supp, supp, supp } ) )
            res += comp( wf.I_b() + I ) * wf.weight(I);

        }, result, field );


    return result;
  }

}

namespace field {
  template < typename T, int DField, int DGrid, typename LocType, typename ShapeF >
  void deposit ( Field<T,DField,DGrid>& field,
                 apt::Vec<T, DField> variable,
                 const LocType& q_abs,
                 const ShapeF& shapef ) {
    const auto& grid = field.mesh().bulk();
    constexpr auto supp = ShapeF::support;
    apt::array<T, DGrid> loc {}; // location at which to interpolate field

    apt::foreach<0,DGrid>
      ( []( auto& l, auto q, const auto& g ){
          l = ( q - g.lower() ) / g.delta();
        }, loc, q_abs, grid );

    apt::foreach<0,DField>
      ( [&] ( const auto& var, auto comp ) { // TODOL comp is proxy, it breaks semantics
          auto wf = WeightFinder( loc, comp.offset, shapef );
          for ( const auto& I : apt::Block( apt::Index<DGrid>{ supp, supp, supp } ) )
            comp( wf.I_b() + I ) += var * wf.weight(I);
        }, variable, field );

    return;
  }

}
