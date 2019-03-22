#include "field/field_shape_interplay.hpp"
#include "apt/block.hpp"
#include "apt/vec.hpp"

// Also see current_deposition.cpp

namespace field::impl {
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

namespace field {

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
