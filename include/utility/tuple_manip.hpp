#ifndef  _TUPLE_MANIP_HPP_
#define  _TUPLE_MANIP_HPP_

#include <cmath>
#include <tuple>
#include <utility>
#include <experimental/type_traits> // for is_detected

namespace tum {
  namespace impl {
    // NOTE only type no varialbe is a common trick to suppress compiler from warning about unused varialbe
    template < typename Unary, typename Tuple, std::size_t... I >
    constexpr auto map_unary( Unary f, Tuple tp, std::index_sequence<I...> ) {
      return std::make_tuple( f(std::get<I>(tp))... );
    }

    template < typename Binary, typename Tuple1, typename Tuple2, std::size_t... I >
    constexpr auto map_binary( Binary f, Tuple1 tp1, Tuple2 tp2, std::index_sequence<I...> ) {
      return std::make_tuple( f(std::get<I>(tp1), std::get<I>(tp2))... );
    }
  }

  template <typename T>
  using is_tuplelike_t = decltype( std::get<std::tuple_size<T>::value - 1>(std::declval<T>()) );

  // in practice, only std::tuple and std::array are tuplelike
  template <typename T>
  constexpr bool is_tuplelike() {
    return std::experimental::is_detected< is_tuplelike_t, T >::value;
  }

  template <typename... T>
  using enable_if_tuple_t = std::enable_if_t< (... && is_tuplelike<T>() ), int>;
}

namespace tum {
  template < typename Unary, class Tuple, enable_if_tuple_t<Tuple> >
  constexpr auto map( Unary&& f, Tuple&& tp ) {
    using Indices = std::make_index_sequence<std::tuple_size<Tuple>::value>;
    // TODOL maybe constexpr lambda can remove the need of unary_impl, but it seems precarious
    return impl::map_unary( std::forward<Unary>(f), std::forward<Tuple>(tp), Indices{} );
  }


  template < typename Binary, typename Tuple1, typename Tuple2, enable_if_tuple_t<Tuple1, Tuple2> >
  constexpr auto map( Binary&&f, Tuple1&& tp1, Tuple2&& tp2 ) {
    static_assert( std::tuple_size<Tuple1>::value == std::tuple_size<Tuple2>::value, "Unmatched size of tuple operands!");
    using Indices = std::make_index_sequence<std::tuple_size<Tuple1>::value>;
    return impl::map_binary( f, tp1, tp2, Indices{} );
  }

}

// TODOL double check += returns current object
// TODOL double check self assignment

template < typename Tuple1, typename Tuple2, tum::enable_if_tuple_t<Tuple1,Tuple2> >
constexpr auto operator+= ( Tuple1&& tp1, Tuple2&& tp2 ) {
  return tum::map( [](auto& lhs, auto rhs) { return lhs += rhs; },
                   std::forward<Tuple1>(tp1), std::forward<Tuple2>(tp2) );
}

template < typename Tuple1, typename Tuple2, tum::enable_if_tuple_t<Tuple1,Tuple2> >
constexpr auto operator-= ( Tuple1&& tp1, Tuple2&& tp2 ) {
  return tum::map( [](auto& lhs, auto rhs) { return lhs -= rhs; },
                   std::forward<Tuple1>(tp1), std::forward<Tuple2>(tp2) );
}

template < typename Tuple, typename T, tum::enable_if_tuple_t<Tuple> >
constexpr auto operator*= ( Tuple&& tp, T a ) {
  return tum::map( [ rhs=a ](auto& lhs) { return lhs *= rhs; },
                   std::forward<Tuple>(tp) );
}

template < typename Tuple, typename T, tum::enable_if_tuple_t<Tuple> >
constexpr auto operator/= ( Tuple&& tp, T a ) {
  return tum::map( [ rhs=a ](auto& lhs) { return lhs /= rhs; },
                   std::forward<Tuple>(tp) );
}

template < typename Tuple1, typename Tuple2, tum::enable_if_tuple_t<Tuple1,Tuple2> >
constexpr auto operator+ ( Tuple1&& tp1, Tuple2&& tp2 ) {
  return tum::map( [](auto lhs, auto rhs) { return lhs + rhs; },
                   std::forward<Tuple1>(tp1), std::forward<Tuple2>(tp2) );
}

template < typename Tuple1, typename Tuple2, tum::enable_if_tuple_t<Tuple1,Tuple2> >
constexpr auto operator- ( Tuple1&& tp1, Tuple2&& tp2 ) {
  return tum::map( [](auto lhs, auto rhs) { return lhs - rhs; },
                   std::forward<Tuple1>(tp1), std::forward<Tuple2>(tp2) );
}

template < typename Tuple, typename T, tum::enable_if_tuple_t<Tuple> >
constexpr auto operator* ( Tuple&& tp, T a ) {
  return tum::map( [ rhs=a ](auto lhs) { return lhs * rhs; },
                   std::forward<Tuple>(tp) );
}

template < typename Tuple, typename T, tum::enable_if_tuple_t<Tuple> >
constexpr auto operator* ( T a, Tuple&& tp ) {
  return std::forward<Tuple>(tp) * a;
}

template < typename Tuple, typename T, tum::enable_if_tuple_t<Tuple> >
constexpr auto operator/ ( Tuple&& tp, T a ) {
  return tum::map( [ rhs=a ](auto lhs) { return lhs / rhs; },
                   std::forward<Tuple>(tp) );
}

namespace tum {
  template < typename Tuple1, typename Tuple2, enable_if_tuple_t<Tuple1,Tuple2> >
  constexpr auto inner_product ( Tuple1&& tp1, Tuple2&& tp2 ) {
    return std::apply( []( auto a, auto b ){return a+b;}, std::forward<Tuple1>(tp1) * std::forward<Tuple2>(tp2) );
  }

  template < typename Tuple, enable_if_tuple_t<Tuple> >
  constexpr auto abs_sq ( Tuple&& tp ) {
    return inner_product( std::forward<Tuple>(tp), std::forward<Tuple>(tp) );
  }

  template < typename Tuple, enable_if_tuple_t<Tuple> >
  constexpr auto abs ( Tuple&& tp ) {
    return std::sqrt( abs_sq( std::forward<Tuple>(tp) ) );
  }
}

#endif
