#include "particle_module/depositer.hpp"

namespace esirkepov :: impl {
  template < std::size_t DGrid, sf::shape S, typename T >
  class ShapeRangeInterator {
  private:
    static constexpr sf::ShapeFunction<S,T> shape_f;

    int _I = 0;

    std::array<T, DGrid - 1> _wgt;
    const Vec<int, DGrid> _I_b;
    const Vec<T, DGrid> _sep_b;

  public:
    using difference_type = int;
    using value_type = void;
    using reference = std::tuple< std::array<int, DGrid>, T >;
    using pointer = void;
    using iterator_category = std::forward_iterator_tag;

    ShapeRangeInterator( int I, const Vec<T,DGrid>& location )
      : _I_b( vec::per_dim::make<DGrid>
              ( []( const auto& loc ) {
                  return int(loc - shape_f.radius) + 1;}, location ) ),
        _sep_b( location - _I_b ) {
    }

    inline bool operator!= ( int I_end ) const {
      return _I != I_end;
    }

    auto& operator++() { ++_I; return *this; }

    // TODO make sure ShapeF can be passed in as constexpr. This may be used to optimize the ever-checking away.
    reference operator*() const {
      constexpr auto supp = sf::support(S);
      int i = _I % supp;
      int j = i / supp;
      int k = i / (supp * supp);
      if ( 0 == i ) {
        wj = shape_f( std::get<1>(sep_b) - j );
        if ( 0 == j ) {
          wk = shape_f( std::get<2>(sep_b) - k );
        }
      }
      return std::make_tuple(_I_b + {i,j,k},
                             shape_f( std::get<0>(sep_b) - i ) * wj * wk );
    }

  };

  template < std::size_t DGrid, sf::shape S, typename T >
  class ShapeRange {
  private:
    Vec<T,DGrid> _loc;
  public:
    ShapeRange( const Vec<T,DGrid>& loc_rel, const Grid<DGrid>& grid, const Vec<T,DGrid>& offset )
      : _loc ( loc_rel - offset - mem::lower(grid) + mem::guard(grid) ) {}

    auto begin() const && {
                           return ShapeRangeInterator<DGrid, S, T>( 0, std::move(_loc) );
    }

  };

  // TODO double check this function simply returns constexpr. We omitted the argument name
  template < std::size_t DGrid, sf::shape S, typename T >
  constexpr std::end( const ShapeRange<DGrid, S, T> & ) noexcept {
    return DGrid * sf::support(S);
  }
}

namespace esirkepov{
  template < std::size_t DGrid, sf::shape S, typename T >
  inline auto make_shape_range( const Vec<T, DGrid>& loc_rel, const Grid<DGrid>& grid, const Vec<T,DGrid>& offset ) {
    return impl::ShapeRange<DGrid, S, T>( loc_rel, grid, offset );
  }


  inline double CalculateW_2DVersion( double sx0, double sx1, double sy0, double sy1 ){
    return ( ( 2 * sx1 + sx0 ) * sy1 + ( sx1 * 2 * sx0 ) * sy0 ) / 6.0;
  }

  inline double CalculateW_3DVersion( double sx0, double sx1,
                                      double sy0, double sy1,
                                      double sz0, double sz1 ){
    return (sx1 - sx0) * CalculateW_2DVersion(sy0, sy1, sz0, sz1);
  }
}

namespace particle {
  namespace esir = esirkepov;
  template < typename Tvt, std::size_t DPtc, std::size_t DField, std::size_t DGrid,
             sf::shape S, typename Trl = vec::remove_cvref_t<Tvt> >
  void depositWJ( const Particle<Tvt,DPtc>& ptc, Field<Trl,DField>& WJ  ) {
    // NOTE static_assert(WJ is the correct stagger)

    const auto& p_tmp = ptc.p;
    Wfactor[2] = p_tmp.z / std::sqrt( 1.0 + p_tmp.dot(p_tmp) );
    std::array<double,3> W = { 0.0, 0.0, 0.0 };
    for ( auto[ I, W ] : esir::make_shape_range() ) {
      WJ[0] += W[0];
      WJ[1] += W[1];
      if constexpr ( DGrid == 2 ) {
          // FIXME fix the following
          // Calling deposition after pusher.calculateDisplacement implies
          // that p_tmp, which is at n+0.5 time step, is based with respect
          // to x^(n+1). However, the x used here is actually x^(n). One way
          // to improve this is obviously calling updatePos before deposition
          // and accordingly change expressions for calculating shapefunctions.
          // FIXME: But, where is J based? Does one really need rebasing momentum?
          const auto& p_tmp = ptc.p;
          WJ[2] += W[2] * p_tmp.z / std::sqrt( 1.0 + p_tmp.dot(p_tmp) )
        } else if ( DGrid == 3 ) {
        WJ[2] += W[2];
      }
    }



    return;
  }
}
