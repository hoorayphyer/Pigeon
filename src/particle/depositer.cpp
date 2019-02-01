#include "particle/depositer.hpp"
#include "apt/numeric.hpp"
#include "kernel/grid.hpp"
#include "kernel/shape.hpp"

#include "field/field.hpp"
#include "particle/particle.hpp"

namespace esirkepov :: impl {
  template < typename T >
  inline T calcW_2D( T sx0, T sx1, T sy0, T sy1 ) noexcept {
    return ( ( 2 * sx1 + sx0 ) * sy1 + ( sx1 * 2 * sx0 ) * sy0 ) / 6.0;
  }

  template < typename T >
  inline T calcW( T sx0, T sx1, T sy0, T sy1, T sz0, T sz1 ) noexcept {
    return (sx1 - sx0) * calcW_2D(sy0, sy1, sz0, sz1);
  }

  template < std::size_t DGrid, knl::shape S, typename T, std::size_t DField,
             class ShapeRange >
  class ShapeRangeInterator {
    static_assert( DField == 3 );
  private:
    int _I = 0;
    const ShapeRange& _sr;

    std::array<int, DGrid> ijk;
    std::array<T, DGrid> s0;
    std::array<T, DGrid> s1;
    std::array<T, DField> W;

  public:
    using difference_type = int;
    using value_type = void;
    using reference = std::tuple< std::array<int, DGrid>, T >;
    using pointer = void;
    using iterator_category = std::forward_iterator_tag;

    ShapeRangeInterator( int I, const ShapeRange& sr ) noexcept : _I(I), _sr(sr) {}

    inline bool operator!= ( const ShapeRangeInterator& other ) const noexcept {
      return _I != other._I;
    }

    auto& operator++() noexcept { ++_I; return *this; }

    // TODO make sure ShapeF can be passed in as constexpr. This may be used to optimize the ever-checking away.
    reference operator*() const noexcept {
      const auto& stride = _sr._stride;

      auto f_calc_s0s1 =
        [&shapef = _sr.shapef] ( auto& s0, auto& s1, const auto& index,
                                 const auto& sep1_b, const auto& dq ) noexcept {
          s0 = shapef( index + sep1_b + dq );
          s1 = shapef( index + sep1_b );
        };

      if constexpr ( DGrid >= 2 ) {
          std::get<0>(ijk) = _I % std::get<0>(stride);
          apt::foreach<0,1>( f_calc_s0s1, s0, s1, ijk, _sr._sep1_b, _sr.dq );

          if ( std::get<0>(ijk) != 0 ) goto CALCW;
          std::get<1>(ijk) = ( _I % std::get<1>(stride) ) / stride[0];
          apt::foreach<1,2>( f_calc_s0s1, s0, s1, ijk, _sr._sep1_b, _sr.dq );

          if constexpr ( DGrid == 3 ) {
              if ( std::get<1>(ijk) != 0 ) goto CALCW;
              std::get<2>(ijk) = _I / stride[1];
              apt::foreach<2,3>( f_calc_s0s1, s0, s1, ijk, _sr._sep1_b, _sr.dq );
            }

        CALCW:

          if constexpr ( DGrid == 2 ) {
              W[0] = calcW( s0[0], s1[0], s0[1], s1[1], 1.0, 1.0 );
              W[1] = calcW( s0[1], s1[1], 1.0, 1.0, s0[0], s1[0] );
              W[2] = calcW( 0.0, 1.0, s0[0], s1[0], s0[1], s1[1] );
            } else {
            W[0] = calcW( s0[0], s1[0], s0[1], s1[1], s0[2], s1[2] );
            W[1] = calcW( s0[1], s1[1], s0[2], s1[2], s0[0], s1[0] );
            W[2] = calcW( s0[2], s1[2], s0[0], s1[0], s0[1], s1[1] );
          }

        } else {
        static_assert( DGrid > 1 && DGrid < 4 );
      }

      return std::forward_as_tuple(_sr._I_b + ijk, W );
    }

  };


  template < std::size_t DGrid, knl::shape S, typename T, std::size_t DField >
  class ShapeRange {
  private:
    const apt::Vec<T,DGrid>& dq;
    const knl::shapef_t<S>& shapef;
    std::array<int, DGrid> _I_b;
    std::array<T, DGrid> _sep1_b;
    std::array<int, DGrid> _stride;

  public:
    friend class ShapeRangeInterator<DGrid, S, T, DField, ShapeRange>;

    ShapeRange( const apt::Vec<T,DGrid>& q1_rel, const apt::Vec<T,DGrid>& dq_rel, const knl::grid_t<DGrid,T>& grid, const knl::shapef_t<S>& shapefunc )
      : dq( dq_rel ), shapef(shapefunc) {
      // NOTE offset is 0.5
      // TODO expr template on grid
      auto&& q1 = q1_rel - 0.5 - mem::lower(grid) + mem::guard(grid);
      apt::foreach<0, DGrid>
        ( []( auto& ind_b, auto& sep1_b, auto& stride, const auto& x1, const auto& dx ) noexcept {
            ind_b = int( std::min(x1, x1-dx) - knl::radius(S) ) + 1;
            sep1_b = ind_b - x1;
            stride = int( std::max(x1, x1-dx) - knl::radius(S) ) + 1 + knl::support(S) - ind_b;
          }, _I_b, _sep1_b, _stride, std::move(q1), dq );

      if constexpr ( DGrid == 2 ) {
        std::get<1>(_stride) *= std::get<0>(_stride);
      } else if ( DGrid == 3 ) {
        std::get<1>(_stride) *= std::get<0>(_stride);
        std::get<2>(_stride) *= std::get<1>(_stride);
      }
      static_assert( DGrid > 1 && DGrid < 4 );

    }

    auto begin() const {
      return ShapeRangeInterator<DGrid, S, T, DField, ShapeRange>( 0, *this );
    }

    auto end() const {
      return ShapeRangeInterator<DGrid, S, T, DField, ShapeRange>( std::get<DGrid-1>(_stride), *this );
    }

  };

}

namespace esirkepov {
  template < std::size_t DGrid, knl::shape S, typename T, std::size_t DField >
  inline auto make_shape_range( const apt::Vec<T, DGrid>& q1_rel, const apt::Vec<T, DGrid>& dq_rel, const knl::grid_t<DGrid,T>& grid ) {
    return impl::ShapeRange<DGrid, S, T, DField>( q1_rel, dq_rel, grid );
  }
}

namespace particle {
  template < knl::shape S, typename T_WJ, typename Tvt,
             std::size_t DPtc, std::size_t DField, std::size_t DGrid
             >
  struct depositWJ_t<S,
                     field::Field<T_WJ,DField,DGrid>,
                     Particle<Tvt,DPtc>,
                     apt::Vec<apt::remove_cvref_t<Tvt>,DPtc>,
                     knl::grid_t<DGrid, apt::remove_cvref_t<Tvt>> > {
    void operator() ( field::Field<T_WJ,DField,DGrid>& WJ, const Particle<Tvt,DPtc>& ptc,
                      const apt::Vec<apt::remove_cvref_t<Tvt>,DPtc>& dq, const knl::grid_t<DGrid, apt::remove_cvref_t<Tvt>>& grid ) {
      namespace esir = esirkepov;
      // NOTE static_assert(WJ is the correct stagger)

      for ( auto[ I, W ] : esir::make_shape_range(ptc.q, dq, grid) ) {
        // TODO optimize indices use. The problem is that I_b + ijk is global, hence when passed to the interface of f( global index ), the same subtraction I_b - anchor will happen many times.
        WJ.c<0>(I) += std::get<0>(W);
        WJ.c<1>(I) += std::get<1>(W);

        // FIXME fix the following in DGrid == 2
        // Calling deposition after pusher.calculateDisplacement implies
        // that p_tmp, which is at n+0.5 time step, is based with respect
        // to x^(n+1). However, the x used here is actually x^(n). One way
        // to improve this is obviously calling updatePos before deposition
        // and accordingly change expressions for calculating shapefunctions.
        // FIXME: But, where is J based? Does one really need rebasing momentum?
        if constexpr ( DGrid == 2 ) {
            WJ.c<2>(I) += std::get<2>(W) * std::get<2>(ptc.p) / std::sqrt( 1.0 + apt::sqabs(ptc.p) );
          } else if ( DGrid == 3 ) {
          WJ.c<2>(I) += std::get<2>(W);
        }
        static_assert( DGrid > 1 && DGrid < 4 );
      }

      return;
    }
  };

}
