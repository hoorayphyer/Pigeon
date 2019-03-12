#ifndef  _APT_TYPE_TRAITS_HPP_
#define  _APT_TYPE_TRAITS_HPP_

namespace apt {
  template < class C >
  inline constexpr int ndim_v = C::NDim;

  template < class C >
  struct element {
    using type = typename C::element_type;
  };

  template < class C >
  using element_t = typename element<C>::type;
}

namespace apt {
  // template < typename... T >
  // using most_precise_t = std::enable_if_t< (... && std::is_arithmetic_v<T>),
  //                                          decltype( (... + (T)0) ) >;
  template < typename... T >
  using most_precise_t = decltype( (... + (T)0) );
}

#endif
