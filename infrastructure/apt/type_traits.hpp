#pragma once

namespace apt {
template <class C>
inline constexpr int ndim_v = C::NDim;

template <class C>
struct element {
  using type = typename C::element_type;
};

template <class C>
using element_t = typename element<C>::type;
}  // namespace apt

namespace apt {
// template < typename... T >
// using most_precise_t = std::enable_if_t< (... && std::is_arithmetic_v<T>),
//                                          decltype( (... + (T)0) ) >;
template <typename... T>
using most_precise_t = decltype((... + (T)0));
}  // namespace apt

namespace apt {
namespace impl {
struct True {
  static constexpr bool value = true;
};
struct False {
  static constexpr bool value = false;
};

template <typename T>
struct is_ref : False {};
template <typename T>
struct is_ref<T&> : True {};
template <typename T>
struct is_ref<T&&> : True {};

template <typename T>
struct is_const : False {};
template <typename T>
struct is_const<const T> : True {};
// TODO std::is_const_v<T> is false if T is a reference
template <typename T>
struct is_const<const T&> : False {};
template <typename T>
struct is_const<const T&&> : False {};

}  // namespace impl
template <typename T>
inline constexpr auto is_ref_v = impl::is_ref<T>::value;
template <typename T>
inline constexpr auto is_const_v = impl::is_const<T>::value;

namespace impl {
template <typename T>
struct remove_ref {
  using type = T;
};
template <typename T>
struct remove_ref<T&> {
  using type = T;
};
template <typename T>
struct remove_ref<T&&> {
  using type = T;
};

// TODO add volatile for completeness
template <typename T>
struct remove_cv {
  using type = T;
};
template <typename T>
struct remove_cv<const T> {
  using type = T;
};
}  // namespace impl
template <typename T>
using remove_ref_t = typename impl::remove_ref<T>::type;
template <typename T>
using remove_cv_t = typename impl::remove_cv<T>::type;

namespace impl {
template <bool Cond, typename True, typename False>
struct cond {
  using type = True;
};

template <typename True, typename False>
struct cond<false, True, False> {
  using type = False;
};
}  // namespace impl
template <bool Cond, typename True, typename False>
using cond_t = typename impl::cond<Cond, True, False>::type;

namespace impl {
template <typename From, typename To>
struct copy_ref {
  using type = cond_t<is_ref_v<From>, To&, To>;
};

template <typename From, typename To>
struct copy_const {
  using type = cond_t<is_const_v<remove_ref_t<From>>, const To, To>;
};
}  // namespace impl

template <typename From, typename To>
using copy_ref_t = typename impl::copy_ref<From, To>::type;

template <typename From, typename To>
using copy_const_t = typename impl::copy_const<From, To>::type;

// TODO add volatile for completeness
template <typename From, typename To>
using copy_cvref_t = copy_const_t<From, copy_ref_t<From, To>>;

// template < typename T, typename U >
// inline constexpr auto is_same_cvref_v = std::is_same_v< apt::copy_cvref_t<T,
// U>, U >;

// NOTE: somehow I don't want to use std::decay. This template will be provided
// in std in C++20 anyway
template <typename T>
using remove_cvref_t = remove_cv_t<remove_ref_t<T>>;
}  // namespace apt
