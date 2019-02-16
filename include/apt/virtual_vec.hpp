#ifndef  _APT_VIRTUAL_VEC_HPP_
#define  _APT_VIRTUAL_VEC_HPP_

#include "apt/vec_slice.hpp"

// NOTE virtual vector, or a vector proxy
namespace apt {
  template < typename T, int N >
  struct vVec : public VecExpression<vVec<T,N>, T> {
  private:
    T& _head;
    vVec<T,N-1>& _tail;

  public:
    using value_type = T;
    static constexpr int size = N;

    template < typename U, typename... R >
    constexpr vVec( U& u, R&... r ) noexcept
      : _head(u), _tail(r...) { static_assert(sizeof...(R) == N - 1); };

    template < int I >
    constexpr T v() const noexcept {
      if constexpr ( I == 0 ) return _head;
      else return _tail.template v<I-1>();
    }

    template < int I >
    constexpr T& v() noexcept {
      if constexpr ( I == 0 ) return _head;
      else return _tail.template v<I-1>();
    }

    // TODO double check this
    constexpr vVec( vVec&& vec ) noexcept
     : _head( std::get<0>(vec) ),
       _tail( VecSlice<1, vVec::size, vVec>(vec) ) {}

    vVec() = delete;
    vVec( const vVec& ) = delete; // because it breaks copy sematics
    template <typename E>
    vVec( const VecExpression<E>& ) = delete;
    template <typename E>
    vVec( VecExpression<E>&& ) = delete;
  };

  template < typename T >
  struct vVec<T,0> {
    template < typename ...Args >
    constexpr vVec( Args&&... ) noexcept {};
  };


}

#endif
