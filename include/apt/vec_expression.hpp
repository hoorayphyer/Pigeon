#ifndef _APT_VEC_EXPRESSION_HPP_
#define _APT_VEC_EXPRESSION_HPP_

#include <experimental/type_traits> // for std::is_detected
#include <type_traits>
#include <utility> // for std::forward

namespace apt {
  template <typename E>
  class VecExpression {
  public:
    static constexpr auto size = E::size;

    // NOTE .template v<I>() is necessary to tell compiler that v is a template. It works like typename
    template < int I >
    constexpr auto v() const noexcept {
      return static_cast<const E&>(*this).template v<I>();
    }
  };
}

namespace std {
  template < int I, typename E >
  constexpr auto get( const apt::VecExpression<E>& vec ) noexcept {
    static_assert( I < E::size );
    return vec.template v<I>();
  }
}

namespace apt {
  template < typename V >
  struct is_lvec { // lvec = lvalued vec
  private:
    template < typename U >
    using copy_assign_t = decltype( std::declval<U&>().template v<0>() = 1.0 );

  public:
    static constexpr bool value = std::experimental::is_detected< copy_assign_t, V >::value;
  };

  template < typename V >
  inline constexpr auto is_lvec_v = is_lvec<V>::value;

  template < typename V >
  using lvec_ref = std::enable_if_t< is_lvec_v<V>, V& >;

}

namespace apt {
  template < int N, typename F, typename... Arg >
  struct VecFromFunction;

  template < int N, typename F, typename Arg >
  struct VecFromFunction<N, F, Arg>
    : public VecExpression<VecFromFunction<N, F, Arg>> {
    static constexpr auto size = N;
    const F& f;
    const Arg& arg;

    template < int I >
    constexpr auto v() const noexcept {
      return f(std::get<I>(arg));
    }

  };

  template < int N, typename F, typename Arg1, typename Arg2 >
  struct VecFromFunction<N, F, Arg1, Arg2>
    : public VecExpression<VecFromFunction<N, F, Arg1, Arg2>> {
    static constexpr auto size = N;
    const F& f;
    const Arg1& arg1;
    const Arg2& arg2;

    template < int I >
    constexpr auto v() const noexcept {
      return f(std::get<I>(arg1), std::get<I>(arg2));
    }

  };


  template < int N, typename F, typename Arg >
  struct vVecFromFunction
    : public VecExpression<VecFromFunction<N, F, Arg>> {
    static constexpr auto size = N;
    const F& f;
    Arg& arg;

    template < int I >
    constexpr auto v() const noexcept {
      return f(std::get<I>(arg));
    }

    template < int I >
    constexpr auto& v() noexcept {
      return f(std::get<I>(arg));
    }
  };

  template < int N, typename F, typename... Arg >
  constexpr auto make_vff( const F& f, const Arg&... arg ) noexcept {
    return VecFromFunction<N, F, Arg...>{f, arg...};
  }

  template < int N, typename F, typename Arg >
  constexpr auto make_vvff( const F& f, Arg& arg ) noexcept {
    return vVecFromFunction<N, F, Arg>{f, arg};
  }

}

namespace apt {
  template < int Begin, int End, typename E >
  struct VecSlice : public VecExpression<VecSlice<Begin, End, E>> {
    static_assert ( 0 <= Begin && Begin < End && End <= E::size  );
    E& e;

    static constexpr auto size = End - Begin;

    template < int I >
    constexpr auto v() const noexcept {
      return e.template v< Begin + I >();
    }

    template < int I >
    constexpr std::enable_if_t< !std::is_const_v<E>, decltype(e.template v< Begin + I >())&>
    v() noexcept {
      return e.template v< Begin + I >();
    }

  };
}


namespace apt {
  template < std::size_t Begin, std::size_t End, typename Func, typename... Vectors >
  constexpr void foreach( const Func& f, Vectors&&... vecs  ) noexcept {
    static_assert( Begin <= End );
    if constexpr ( Begin == End ) return;
    else {
      f( std::get<Begin>(std::forward<Vectors>(vecs))... );
      return foreach<Begin+1, End>( f, std::forward<Vectors>(vecs)... );
    }
  }

}

#endif
