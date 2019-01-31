#ifndef  _APT_VEC_HPP_
#define  _APT_VEC_HPP_

#include <array>
#include <tuple>
#include <type_traits>

namespace apt {
  template < typename T, std::size_t N >
  struct homogeneous_tuple_gen {
    using type = decltype( std::tuple_cat( std::declval<std::tuple<T>>(), std::declval<typename homogeneous_tuple_gen<T, N-1>::type>() ) );
  };

  template < typename T >
  struct homogeneous_tuple_gen <T,0> { using type = std::tuple<>; };

  template < typename T, std::size_t N >
  using homogeneous_tuple = typename homogeneous_tuple_gen<T,N>::type;

  // NOTE tuple<void> is allowed but cannot be used to instantiate. Its use is kinda restricted to querying it's type like std::tuple_element_t<0, tuple<void>>
  // NOTE lambda expression is not allowed in unevaluated context, such as decltype( [](){...} )

  template < typename Trl, std::size_t N >
  struct Vec : public std::array<Trl,N> {
    template < typename U >
    constexpr Vec( std::array<U,N> x ) noexcept
      : std::array<Trl,N>( std::move(x) ) {}

    // conversion from a virtual vector
    constexpr Vec( const Vec<Trl&,N>& other ) noexcept;
  };


  template < typename Trl, std::size_t N >
  struct Vec<Trl&, N> : public homogeneous_tuple<Trl&,N> {
    constexpr Vec( homogeneous_tuple<T&,N> x ) noexcept
      : homogeneous_tuple<T&,N>( std::move(x) ) {}

    // TODO double check the semantics of these
    Vec( const Vec& other ) = default;
    Vec( Vec&& other ) noexcept = default;
  };

  template < typename Trl, std::size_t N >
  constexpr Vec<Trl,N>::Vec( const Vec<Trl&,N> & vec ) noexcept {
    std::get<0>(*this) = std::get<0>(vec);
    if constexpr ( N > 1 ) std::get<1>(*this) = std::get<1>(vec);
    if constexpr ( N > 2 ) std::get<2>(*this) = std::get<2>(vec);
    static_assert( N < 4 );
  }
}

namespace std {
  // this section is to make apt::Vec tuple-like
  // NOTE when a derived class is used as a template argument, there is just no conversion to its base class ever.
  template < std::size_t I, typename T, std::size_t N >
  struct tuple_element < I, apt::Vec<T,N> > {
    using type = T;
  };

  template < typename T, std::size_t N >
  struct tuple_size < apt::Vec<T,N> > {
    static constexpr auto value = N;
  };
}

#include <experimental/array> // for make_array
namespace apt {
  template < class Tuple >
  using element_t = std::tuple_element_t<0, std::remove_reference_t<Tuple> >;

  template < class Tuple >
  inline constexpr auto size_v = std::tuple_size_v< std::remove_reference_t<Tuple> >;

  template < typename... Args >
  constexpr auto make( Args&&... args ) noexcept {
    auto&& tmp = std::experimental::make_array( std::forward<Args>(args)... );
    return Vec< element_t<decltype(tmp)>, size_v<decltype(tmp)> >( std::move(tmp) );
  }

  template < typename... Args >
  constexpr auto tie( Args&... args ) noexcept {
    auto&& tmp = std::tie( args... );
    return Vec< element_t<decltype(tmp)>, size_v<decltype(tmp)> >( std::move(tmp) );
  }


}

#endif
