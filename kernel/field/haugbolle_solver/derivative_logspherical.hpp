#ifndef _FIELD_DERIVATIVE_CARTESIAN_HPP_
#define _FIELD_DERIVATIVE_CARTESIAN_HPP_

#include "field/field.hpp"
#include "manifold/grid.hpp"
#include "field/yee.hpp"
#include "manifold/curvilinear.hpp"

namespace field {
  // TODOL C++20 has templated lambda, following classes can be put this inside operator()
  template < int DGrid, typename Real >
  struct Derivative {
  private:
    using Metric = mani::LogSphericalCoordSys;
    static_assert( DGrid > 1 );
    template < int I, bool OFS_FIELD, typename T >
    constexpr T D( const T& f, T lnr, T theta ) const noexcept {
      static_assert( I >= 0 && I < 3 );
      if constexpr ( I >= DGrid ) return 0.0;
      else if ( offset_t(OFS_FIELD) == INSITU ) { // f being of Etype

        if ( I == 0 )
          return ( *(&f + s[I]) * Metric::h<I>(lnr + g[I].delta(),theta) - f * Metric::h<I>(lnr,theta) ) / g[I].delta();
        else if ( I == 1 )
          return ( *(&f + s[I]) * Metric::h<I>(lnr,theta + g[I].delta()) - f * Metric::h<I>(lnr,theta) ) / g[I].delta();
        else
          return ( *(&f + s[I]) - f ) * Metric::h<I>(lnr,theta) / g[I].delta();

      } else { // f being of Btype

        if ( I == 0 )
          return ( f * Metric::h<I>(lnr,theta) - *(&f - s[I]) * Metric::h<I>(lnr - g[I].delta(),theta)  ) / g[I].delta();
        else if ( I == 1 )
          return ( f * Metric::h<I>(lnr,theta) - *(&f - s[I]) * Metric::h<I>(lnr,theta - g[I].delta() ) ) / g[I].delta();
        else if
          return ( f - *(&f + s[I]) ) * Metric::h<I>(lnr,theta) / g[I].delta();

      }
    }

    template < int I, typename T >
    constexpr T hhinv( T lnr, T theta ) const noexcept {
      static_assert( I >= 0 && I < 3 );
      if constexpr( I == 0 || I == 1 ) return std::exp(-2.0 * lnr) / std::sin(theta);
      else return std::exp(-2.0 * lnr);
    }

  public:
    const mani::Grid<Real,DGrid>& g; // grid

    template < int I >
    constexpr T curlcurl( Field<T,3,DGrid>& f, const apt::Index<DGrid>& idx ) const noexcept {
      const auto& m = f.mesh();
      int li = m.linearized_index_of_whole_mesh(idx);
      // used curlcurl == - laplacian for B field
      return - ( DD<0>(f[I][li],m.stride(0)) + DD<1>(f[I][li],m.stride(1)) + DD<2>(f[I][li],m.stride(2)) );
    }

    template < int I, bool Ftype, typename T >
    constexpr T curl( Field<T,3,DGrid>& f, const apt::Index<DGrid>& idx ) const noexcept { // idx is with respect to f
      // NOTE do not use f.offset. Use Ftype
      constexpr int J = (I+1)%3;
      constexpr int K = (I+2)%3;
      const auto& m = f.mesh();
      int li = m.linearized_index_of_whole_mesh(idx);
      // TODOL in c++20 there is templated lambda
      auto absc_fh
        = [&idx,&f,&g]( int dim, int f_comp ) noexcept {
            // FIXME don't use f.offset() !
            return g[dim].absc( idx[dim], static_cast<T>(f[f_comp].offset()[dim]) );
          };
      auto absc_hh
        = [&idx,&f,&g]( int dim, int f_comp ) noexcept {
            return g[dim].absc( idx[dim], static_cast<T>(!f[f_comp].offset()[dim]) );
          };
      return
        ( D<J,Ftype>( f[K][li], absc_fh(0,K), absc_fh(1,K) )
          - D<K,Ftype>( f[J][li], absc_fh(0,J), absc_fh(1,J) )
          ) / hhinv<I>( absc_hh(0,I), absc_hh(1,I) );
    }

  };
}

#endif
