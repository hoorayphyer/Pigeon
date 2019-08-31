#ifndef _FIELD_DERIVATIVE_CARTESIAN_HPP_
#define _FIELD_DERIVATIVE_CARTESIAN_HPP_

#include "field/field.hpp"
#include "manifold/grid.hpp"
#include "field/yee.hpp"

namespace field {
  // TODOL C++20 has templated lambda, following classes can be put this inside operator()
  template < int DGrid, typename T >
  struct Derivative {
  private:
    constexpr T diff_plus( const T& f, int stride ) const noexcept {
      return *(&f + stride) - f;
    }

    constexpr T diff_minus( const T& f, int stride ) const noexcept {
      return  f - *(&f - stride);
    }

    constexpr T diff2( const T& f, int stride ) const noexcept {
      return *(&f + stride) + *(&f - stride) - 2.0 * f;
    }

    template < int I, offset_t OFS_FIELD >
    constexpr T D( const T& f, int stride ) const noexcept {
      static_assert( I >= 0 && I < 3 );
      if constexpr ( I >= DGrid ) return 0.0;
      else if ( offset_t(OFS_FIELD) == INSITU ) { // f being of Etype
        return diff_plus(f,stride) / g[I].delta();
      } else { // f being of Btype
        return diff_minus(f,stride) / g[I].delta();
      }
    }

    template < int I >
    constexpr T DD( const T& f, int stride ) const noexcept {
      if constexpr ( I >= DGrid ) return 0.0;
      else return diff2( f, stride ) / (g[I].delta() * g[I].delta());
    }

  public:
    const mani::Grid<T,DGrid>& g; // grid

    template < int I >
    constexpr T curlcurl( Field<T,3,DGrid>& f, const apt::Index<DGrid>& idx ) const noexcept {
      const auto& m = f.mesh();
      int li = m.linearized_index_of_whole_mesh(idx);
      // used curlcurl == - laplacian for B field
      return - ( DD<0>(f[I][li],m.stride(0)) + DD<1>(f[I][li],m.stride(1)) + DD<2>(f[I][li],m.stride(2)) );
    }

    template < int I, offset_t Ftype >
    constexpr T curl( Field<T,3,DGrid>& f, const apt::Index<DGrid>& idx ) const noexcept { // idx is with respect to f
      // NOTE do not use f.offset. Use Ftype
      constexpr int J = (I+1)%3;
      constexpr int K = (I+2)%3;
      const auto& m = f.mesh();
      int li = m.linearized_index_of_whole_mesh(idx);
      // TODOL in c++20 there is templated lambda
      return D<J,Ftype>( f[K][li], m.stride(J) ) - D<K,Ftype>( f[J][li], m.stride(K) );
    }

  };
}

#endif
