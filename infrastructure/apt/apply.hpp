#pragma once

#include <utility>  // for std::forward

namespace apt {
namespace impl {
template <typename Func, typename Arg, int... I>
constexpr decltype(auto) apply_impl(const Func& f, Arg&& arg,
                                    std::integer_sequence<I...>) noexcept {
  // TODO add bounds check
  // static_assert( sizeof...(I) <= Arg::NDim );
  return f(std::forward<Arg>(arg)[I]...);
}
}  // namespace impl

template <typename Func, typename Arg, int D = Arg::NDim>
constexpr decltype(auto) apply(const Func& f, Arg&& arg) noexcept {
  return impl::apply_impl(f, std::forward<Arg>(arg),
                          std::make_integer_sequence<int, D>{});
}
}  // namespace apt
