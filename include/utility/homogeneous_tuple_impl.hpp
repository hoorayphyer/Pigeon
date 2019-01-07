#ifndef  _HOMOGENEOUS_TUPLE_IMPL_HPP_
#define  _HOMOGENEOUS_TUPLE_IMPL_HPP_

#include <tuple>
#include <type_traits>

namespace hmt :: impl {
  template < typename T, std::size_t N >
  struct tuple_concat {
    // TODO check std::declval type
    using type = decltype( std::tuple_cat( std::declval<std::tuple<T&&>>(), std::declval<typename tuple_concat<T&&, N-1>::type>() ) );
  };

  template < typename T >
  struct tuple_concat<T&&,0> { using type = std::tuple<>; };

}

namespace hmt :: impl {

  template < typename T, std::size_t N,
             class Base = typename tuple_concat<T&&,N>::type,
             class = std::enable_if_t< (N > 0), int > >
  struct homotuple : public Base {
  private:
    explicit homotuple( Base tp ) noexcept : Base( std::move(tp) ) {}
  public:
    template<typename ... Types>
    friend auto homogeneitize( std::tuple<Types...>&& ) noexcept;
  };

  template < typename... Types >
  constexpr bool is_homogeneous() noexcept {
    static_assert( sizeof...(Types) > 0 );
    using T_1st = typename std::tuple_element< 0, std::tuple<Types...> >::type;
    return (... && std::is_same_v< T_1st, Types >);
  }

  template < typename... Types >
  inline auto homogeneitize( std::tuple<Types...>&& tp ) noexcept {
    static_assert( is_homogeneous<Types...>(), " all argument types must be the same");

    using T_1st = typename std::tuple_element< 0, std::tuple<Types...> >::type;
    constexpr auto size = sizeof...(Types);
    return homotuple<T_1st, size>( std::move(tp) );
  }

  template < typename... Args >
  inline auto make_homotuple( Args&&... args ) {
    return homogeneitize( std::make_tuple( std::forward<Args>(args)...) );
  }

  template < typename... Args >
  inline auto tie_to_homotuple( Args&&... args ) noexcept {
    return homogeneitize( std::tie( std::forward<Args>(args)...) );
  }

  template < typename... Args >
  inline auto forward_as_homotuple( Args&&... args ) noexcept {
    return homogeneitize( std::forward_as_tuple( std::forward<Args>(args)...) );
  }
}

#endif
