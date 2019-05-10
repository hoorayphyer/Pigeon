#include "msh/mesh_shape_interplay.hpp"
#include "apt/block.hpp"

namespace msh::impl {
  template < typename T, int DGrid, int Dloc, typename ShapeF >
  struct WeightFinder {
  private:
    static_assert(Dloc >= DGrid);
    apt::Index<DGrid> _I_b {};
    apt::array<T,DGrid> _sep_b {};
    apt::array<T,DGrid> _wgt {};
    const ShapeF& _shapef;
  public:
    // NOTE loc is the relative index
    constexpr WeightFinder( const apt::array<T, Dloc>& loc,
                            const apt::array< field::offset_t, DGrid >& offset,
                            const ShapeF& shapef ) noexcept
      : _shapef( shapef ) {

      apt::foreach<0,DGrid>
        ( []( auto& i_b, auto& s_b, auto l, const auto& ofs ) {
            l -= ofs; // now l is the native grid index

            static_assert( ShapeF::support() > 1); // TODOL

            auto int_flr =
              []( auto q ) noexcept {
                // Since min(q0_std) = 0.0 by design, min(q_nat) = - r - offset, so q_nat + ( support ) / 2.0 + offset >= 0. We will simply use (support + 1) as the shift
                constexpr auto shift = 1 + ShapeF::support();
                return int(q+shift) - shift; };

            i_b = int_flr(l - ShapeF::support() / 2.0) + 1; // i_b is native as well as standard
            s_b = i_b - l;
          }, _I_b, _sep_b, loc, offset );
    }

    constexpr const auto& I_b() const noexcept { return _I_b; }

    constexpr T weight ( const apt::Index<DGrid>& I ) noexcept{
      T result = 1.0; // NOTE 1.0, not 0.0

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
      update_wgt(I);
      apt::foreach<0,DGrid>
        ( [&result]( auto w ){ result *= w; }, _wgt );
      return result;
    };
  };
}

namespace msh {

  template < int DGrid, int Supp >
  constexpr apt::Index<DGrid> ShapeExtent
  ( [](){
      apt::Index<DGrid> res;
      for ( int i = 0; i < DGrid; ++i ) res[i] = Supp;
      return res;}() );

  template < typename T, int DGrid, int Dq, typename ShapeF >
  T interpolate ( const typename field::Component<T,DGrid>& fcomp,
                  const apt::array<T,Dq>& q_std,
                  const ShapeF& shapef ) noexcept {
    T res{};

    auto wf = impl::WeightFinder( q_std, fcomp.offset(), shapef );
    for ( const auto& I : apt::Block(ShapeExtent<DGrid, ShapeF::support()>) )
      res += fcomp( wf.I_b() + I ) * wf.weight(I);

    return res;
  }

  template < typename T, int DField, int DGrid, int Dq, typename ShapeF >
  apt::Vec<T, DField> interpolate ( const field::Field<T,DField,DGrid>& field,
                                    const apt::array<T,Dq>& q_std,
                                    const ShapeF& shapef ) noexcept {
    apt::Vec<T, DField> result;
    constexpr auto supp = ShapeF::support();

    apt::foreach<0,DField>
      ( [&] ( auto& res, const auto& comp ) {
          res = interpolate(comp, q_std, shapef );
        }, result, field );

    return result;
  }

}

namespace msh {
  template < typename T, int DField, int DGrid, int Dq, typename ShapeF >
  void deposit ( field::Field<T,DField,DGrid>& field,
                 apt::array<T, DField> variable,
                 const apt::array<T, Dq>& q_std,
                 const ShapeF& shapef ) noexcept {
    constexpr auto supp = ShapeF::support();

    constexpr apt::Index<DGrid> ext
      ( [supp](){
          apt::Index<DGrid> res;
          for ( int i = 0; i < DGrid; ++i ) res[i] = supp;
          return res;}() );

    apt::foreach<0,DField>
      ( [&] ( const auto& var, auto comp ) { // TODOL comp is proxy, it breaks semantics
          auto wf = impl::WeightFinder( q_std, comp.offset(), shapef );
          for ( const auto& I : apt::Block( ext ) )
            comp( wf.I_b() + I ) += var * wf.weight(I);
        }, variable, field );

    return;
  }

}
