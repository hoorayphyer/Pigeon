#ifndef APT_VEC_FROM_FUNCTION_HPP_
#define APT_VEC_FROM_FUNCTION_HPP_

#include "apt/vec_expression.hpp"
#include <functional>

namespace apt {
  template < int N, typename F, typename... Args >
  struct VecFromFunction
    : public VecExpression<VecFromFunction<N, F, Args...>,
                           std::invoke_result_t<F, Args&&...> > {
  private:
    const F& _f;
    std::tuple<Args&&...> _args;

  public:
    using value_type = std::invoke_result_t<F, Args&&...>;
    static constexpr int size = N;

    constexpr VecFromFunction( const F& f, Args&&... args ) noexcept
      : _f(f), _args( std::forward<Args>(args)... ) {}

    template < int I >
    constexpr value_type v() noexcept {
      return std::invoke( _f, _args );
    }
  };

  template < int N, typename F, typename... Args >
  constexpr auto make_vff( const F& f, Args&&... args ) noexcept {
    return VecFromFunction<N, F, Args...>(f, std::forward<Args>(args)...);
  }

}

#endif
