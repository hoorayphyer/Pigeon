#ifndef APT_VEC_FROM_FUNCTION_HPP_
#define APT_VEC_FROM_FUNCTION_HPP_

#include "apt/vec_expression.hpp"
#include <functional> // std::invoke

namespace apt :: impl {
  template < typename ComponentWiseFunc, typename... VecArgs >
  struct invoke_result {
    using type = std::invoke_result_t<
      ComponentWiseFunc,
      decltype(std::get<0>(std::declval<VecArgs&&>()))...
      >;
  };

  template < typename ComponentWiseFunc, typename... VecArgs >
  using invoke_result_t = typename invoke_result<ComponentWiseFunc, VecArgs...>::type;
}

namespace apt {
  // TODO check this
  template < int N, typename ComponentWiseFunc, typename... VecArgs >
  struct VecFromFunction
    : public VecExpression<VecFromFunction<N, ComponentWiseFunc, VecArgs...>,
                           impl::invoke_result_t<ComponentWiseFunc, VecArgs...> > {
  private:
    const ComponentWiseFunc& _f;
    std::tuple<VecArgs&&...> _vec_args;

    template < int Comp, std::size_t... VecNo_I >
    constexpr decltype(auto) v( std::index_sequence<VecNo_I...> ) noexcept {
      return std::invoke( _f, std::get<Comp>(std::get<VecNo_I>(_vec_args))... );
    }

  public:
    // TODO check if value_type need to be decayed
    using value_type = impl::invoke_result_t<ComponentWiseFunc, VecArgs...>;
    static constexpr int size = N;

    constexpr VecFromFunction( const ComponentWiseFunc& f, VecArgs&&... args ) noexcept
      : _f(f), _vec_args( std::forward<VecArgs>(args)... ) {}

    template < int I >
    constexpr decltype(auto) v() noexcept {
      return v<I>( std::index_sequence_for<VecArgs...>{} );
    }
  };

  template < int N, typename ComponentWiseFunc, typename... Args >
  constexpr auto make_vff( const ComponentWiseFunc& f, Args&&... args ) noexcept {
    return VecFromFunction<N, ComponentWiseFunc, Args...>(f, std::forward<Args>(args)...);
  }

}

#endif
