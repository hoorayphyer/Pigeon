#ifndef  _APT_TYPE_TRAITS_HPP_
#define  _APT_TYPE_TRAITS_HPP_

#include "apt/vec.hpp"

namespace apt {
  namespace impl {
    template < typename T >
    struct is_vec { static constexpr bool value = false; };

    template < typename T, std::size_t N >
    struct is_vec< Vec<T,N> > { static constexpr bool value = true; };
    template < typename T, std::size_t N >
    struct is_vec< Vec<T,N>& > { static constexpr bool value = true; };
    template < typename T, std::size_t N >
    struct is_vec< const Vec<T,N>& > { static constexpr bool value = true; };
    template < typename T, std::size_t N >
    struct is_vec< Vec<T,N>&& > { static constexpr bool value = true; };
    template < typename T, std::size_t N >
    struct is_vec< const Vec<T,N>&& > { static constexpr bool value = true; };
  }

  template < typename T >
  inline constexpr bool is_vec_v = impl::is_vec<T>::value;
}

namespace apt {
  namespace impl {
    template < typename From, typename To >
    struct copy_ref {
      using type = std::conditional_t< std::is_reference_v<From>, To&, To >;
    };

    // TODO std::is_const_v<T> is false if T is a reference
    template < typename From, typename To >
    struct copy_const {
      using type = std::conditional_t< std::is_const_v< std::remove_reference_t<From> >, const To, To >;
    };

  }

  template < typename From, typename To >
  using copy_ref_t = typename copy_ref<From,To>::type;

  template < typename From, typename To >
  using copy_const_t = typename copy_const<From,To>::type;

  // TODO add volatile for completeness
  template < typename From, typename To >
  using copy_cvref_t = typename copy_const_t<From, copy_ref_t<From, To>>;

  template < typename T, typename U >
  using is_same_cvref_v = std::is_same_v< apt::copy_cvref_t<T, U>, U>;

  // NOTE: somehow I don't want to use std::decay. This template will be provided in std in C++20 anyway
  template < typename T >
  using remove_cvref_t = std::remove_cv_t<std::remove_reference_t<T>>;
}

#endif
