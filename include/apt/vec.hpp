#ifndef  _APT_VEC_HPP_
#define  _APT_VEC_HPP_

#include "apt/vec_expression.hpp"
#include <array>

namespace apt {
  template < typename T, int N >
  struct Vec : public VecExpression<Vec<T,N>> {
  private:
    std::array<T,N> _v;

  public:
    static constexpr auto size = N;

    template < typename U >
    constexpr Vec( U x[N] ) noexcept {
      foreach<0,N>( [](auto& a, const auto& b){ a = b;}, _v, x);
    }

    template < typename E >
    constexpr Vec( const VecExpression<E>& vec ) noexcept {
      foreach<0,N>( [](auto& a, const auto& b){ a = b;}, _v, vec);
    }

    template < typename E >
    constexpr Vec( VecExpression<E>&& vec ) noexcept {
      foreach<0,N>( [](auto& a, auto& b){ std::swap(a,b); }, _v, vec );
    }

    template < int I >
    constexpr T& v() noexcept {
      return std::get<I>(_v);
    }

  };
}


namespace apt {
  // NOTE virtual vector, or a vector proxy
  template < typename T, int N >
  struct vVec : public VecExpression<vVec<T,N>> {
  private:
    T& _head;
    vVec<T,N-1>& _tail;

  public:
    static constexpr auto size = N;

    template < typename U, typename... R >
    constexpr vVec( U& u, R&... r ) noexcept
      : _head(u), _tail(r...) { static_assert(sizeof...(R) == N - 1); };

    template < int I >
    constexpr double v() const noexcept {
      static_assert( I < N );
      if constexpr ( I == 0 ) return _head;
      else return _tail.template v<I-1>();
    }

    template < int I >
    constexpr T& v() noexcept {
      static_assert( I < N );
      if constexpr ( I == 0 ) return _head;
      else return _tail.template v<I-1>();
    }

    // TODO double check this
    vVec( vVec&& vec ) noexcept
     : _head( std::get<0>(vec) ),
       _tail( VecSlice<1, vVec::size, vVec>(vec) ) {}

    vVec() = delete;
    vVec( const vVec& ) = delete; // because it breaks copy sematics
    template <typename E>
    vVec( const VecExpression<E>& ) = delete;
    template <typename E>
    vVec( VecExpression<E>&& ) = delete;
  };

  // TODO fix the ctors
  template < typename T >
  struct vVec<T,1> : public VecExpression<vVec<T,1>> {
  private:
    T& _head;

  public:
    template < typename U >
    constexpr vVec( U& u ) noexcept : _head(u) {}

    template <typename E>
    constexpr vVec( VecExpression<E>& vec ) noexcept
      : _head(std::get<0>(vec)) {}

    template < int I >
    constexpr T v() const noexcept {
      static_assert( I == 0 );
      return _head;
    }

    template < int I >
    constexpr T& v() noexcept {
      static_assert( I == 0 );
      return _head;
    }
  };


}


// namespace std {
//   // this section is to make apt::Vec tuple-like
//   // NOTE when a derived class is used as a template argument, there is just no conversion to its base class ever.
//   template < std::size_t I, typename T, std::size_t N >
//   struct tuple_element < I, apt::Vec<T,N> > {
//     using type = T;
//   };

//   template < typename T, std::size_t N >
//   struct tuple_size < apt::Vec<T,N> > {
//     static constexpr auto value = N;
//   };
// }


#endif
