#ifndef _FIELD_DERIVATIVE_CARTESIAN_HPP_
#define _FIELD_DERIVATIVE_CARTESIAN_HPP_

#include "field/field.hpp"
#include "manifold/grid.hpp"
#include "field/yee.hpp"
#include "manifold/curvilinear_impl.hpp"

namespace field {
  // TODOL C++20 has templated lambda, following classes can be put this inside operator()
  template < int DGrid, typename T >
  struct Derivative {
  private:
    using Metric = mani::LogSphericalCoordSys;
    static_assert( DGrid > 1 );
    template < int I, offset_t OFS_FIELD >
    constexpr T D( const T& f, T lnr, T theta ) const noexcept {
      static_assert( I >= 0 && I < 3 );
      if constexpr ( I >= DGrid ) return 0.0;
      else if ( offset_t(OFS_FIELD) == INSITU ) { // f being of Etype

        if ( I == 0 )
          return ( *(&f + s[I]) * Metric::template h<I>(lnr + g[I].delta(),theta) - f * Metric::template h<I>(lnr,theta) ) / g[I].delta();
        else if ( I == 1 )
          return ( *(&f + s[I]) * Metric::template h<I>(lnr,theta + g[I].delta()) - f * Metric::template h<I>(lnr,theta) ) / g[I].delta();
        else
          return ( *(&f + s[I]) - f ) * Metric::template h<I>(lnr,theta) / g[I].delta();

      } else { // f being of Btype

        if ( I == 0 )
          return ( f * Metric::template h<I>(lnr,theta) - *(&f - s[I]) * Metric::template h<I>(lnr - g[I].delta(),theta)  ) / g[I].delta();
        else if ( I == 1 )
          return ( f * Metric::template h<I>(lnr,theta) - *(&f - s[I]) * Metric::template h<I>(lnr,theta - g[I].delta() ) ) / g[I].delta();
        else
          return ( f - *(&f - s[I]) ) * Metric::template h<I>(lnr,theta) / g[I].delta();

      }
    }

    template < int I >
    constexpr T hhinv( T lnr, T theta ) const noexcept {
      static_assert( I >= 0 && I < 3 );
      if constexpr( I == 0 || I == 1 ) return std::exp(-2.0 * lnr) / std::sin(theta);
      else return std::exp(-2.0 * lnr);
    }

  public:
    const mani::Grid<T,DGrid>& g; // grid
    const apt::Index<DGrid> s; // stride

    template < int I, offset_t Ftype >
    constexpr T curl( Field<T,3,DGrid>& f, const apt::Index<DGrid>& idx ) const noexcept { // idx is with respect to f
      // NOTE do not use f.offset. Use Ftype
      constexpr int J = (I+1)%3;
      constexpr int K = (I+2)%3;
      int li = f.mesh().linearized_index_of_whole_mesh(idx);
      // TODOL in c++20 there is templated lambda
      auto absc_fh
        = [&idx,&g=this->g]( int dim, int f_comp ) noexcept {
            // needs f's offset, which is Ftype when dim == f_comp
            return g[dim].absc( idx[dim], 0.5 * yee::ofs_get<Ftype>(f_comp,dim) );
          };
      auto absc_hh
        = [&idx,&g=this->g]( int dim, int f_comp ) noexcept {
            // needs f's anti-offset, which is !Ftype when dim == f_comp
            return g[dim].absc( idx[dim], 0.5 * yee::ofs_get<!Ftype>(f_comp,dim) );
          };
      return
        ( D<J,Ftype>( f[K][li], absc_fh(0,K), absc_fh(1,K) )
          - D<K,Ftype>( f[J][li], absc_fh(0,J), absc_fh(1,J) )
          ) / hhinv<I>( absc_hh(0,I), absc_hh(1,I) );
    }

  };
}

#endif
