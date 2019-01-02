#ifndef  _TUPLE_MANIP_HPP_
#define  _TUPLE_MANIP_HPP_

#include <cmath>
#include <tuple>
#include <utility>
#include <experimental/type_traits> // for is_detected

namespace tum {

  template <typename T>
  using is_tuplelike_t = decltype( std::get<0>(std::declval<T>()) );

  // in practice, only std::tuple and std::array are tuplelike
  template <typename T>
  constexpr bool is_tuplelike() {
    return std::experimental::is_detected< is_tuplelike_t, T >::value;
  }

  template <typename... T>
  using enable_if_tuple_t = std::enable_if_t< (... && is_tuplelike<T>() ), int>;
  // using enable_if_tuple_t = int;

  template <typename Tuple, class = enable_if_tuple_t<Tuple>>
  constexpr auto tuple_size() {
    return std::tuple_size<typename std::remove_reference<Tuple>::type>::value;
  }


  template < typename Unary, typename Tuple, std::size_t... I >
  constexpr auto map_impl_unary( Unary f, Tuple tp, std::index_sequence<I...> ) {
    return std::make_tuple( f(std::get<I>(tp))... );
  }

  template < typename Func, class... Tuple, std::size_t... I >
  constexpr auto map_impl( Func&& f, std::tuple<Tuple... > tps, std::index_sequence<I...> ) {
    return std::make_tuple( std::apply(f, map_impl_unary( [](auto&& x){return std::get<I>(x);}, tps, std::make_index_sequence<sizeof...(Tuple)>{} ) )... );

  }

  template < typename Func, class Tuple, class... RestTuple >
  constexpr auto map( Func&& f, Tuple&& tp, RestTuple&&... rest_tp ) {
    using Indices = std::make_index_sequence<tuple_size<Tuple>()>;
    if constexpr( sizeof...(rest_tp) == 0 ) {
      return  map_impl_unary( std::forward<Func>(f), std::forward<Tuple>(tp), Indices{} );
    } else {
      // TODO choose the minimum size insteal of forcing all sizes to be equal. Or requiring size of RestTuple be >= size of Tuple. Needed in inner product
      static_assert( (... && (tuple_size<RestTuple>() >= tuple_size<Tuple>()) ), "Unmatched size of tuple operands!");
      return map_impl( std::forward<Func>(f), std::forward_as_tuple(tp, rest_tp...), Indices{} );
    }
  }
}

// TODOL double check += returns current object
// TODOL double check self assignment

template < typename Tuple1, typename Tuple2, class = tum::enable_if_tuple_t<Tuple1,Tuple2> >
constexpr auto operator+= ( Tuple1&& tp1, Tuple2&& tp2 ) {
  return tum::map( [](auto& lhs, auto rhs) { return lhs += rhs; },
                   std::forward<Tuple1>(tp1), std::forward<Tuple2>(tp2) );
}

template < typename Tuple1, typename Tuple2, class = tum::enable_if_tuple_t<Tuple1,Tuple2> >
constexpr auto operator-= ( Tuple1&& tp1, Tuple2&& tp2 ) {
  return tum::map( [](auto& lhs, auto rhs) { return lhs -= rhs; },
                   std::forward<Tuple1>(tp1), std::forward<Tuple2>(tp2) );
}

template < typename Tuple, typename T, class = tum::enable_if_tuple_t<Tuple> >
constexpr auto operator*= ( Tuple&& tp, T a ) {
  return tum::map( [ rhs=a ](auto& lhs) { return lhs *= rhs; },
                   std::forward<Tuple>(tp) );
}

template < typename Tuple, typename T, class = tum::enable_if_tuple_t<Tuple> >
constexpr auto operator/= ( Tuple&& tp, T a ) {
  return tum::map( [ rhs=a ](auto& lhs) { return lhs /= rhs; },
                   std::forward<Tuple>(tp) );
}

template < typename Tuple1, typename Tuple2, class = tum::enable_if_tuple_t<Tuple1,Tuple2> >
constexpr auto operator+ ( const Tuple1& tp1, const Tuple2& tp2 ) {
  return tum::map( [](auto lhs, auto rhs) { return lhs + rhs; }, tp1, tp2 );
}

template < typename Tuple1, typename Tuple2, class = tum::enable_if_tuple_t<Tuple1,Tuple2> >
constexpr auto operator- ( Tuple1&& tp1, Tuple2&& tp2 ) {
  return tum::map( [](auto lhs, auto rhs) { return lhs - rhs; },
                   std::forward<Tuple1>(tp1), std::forward<Tuple2>(tp2) );
}

template < typename Tuple, typename T, class = tum::enable_if_tuple_t<Tuple> >
constexpr auto operator* ( Tuple&& tp, T a ) {
  return tum::map( [ rhs=a ](auto lhs) { return lhs * rhs; },
                   std::forward<Tuple>(tp) );
}

template < typename Tuple, typename T, class = tum::enable_if_tuple_t<Tuple> >
constexpr auto operator* ( T a, Tuple&& tp ) {
  return std::forward<Tuple>(tp) * a;
}

template < typename Tuple, typename T, class = tum::enable_if_tuple_t<Tuple> >
constexpr auto operator/ ( Tuple&& tp, T a ) {
  return tum::map( [ rhs=a ](auto lhs) { return lhs / rhs; },
                   std::forward<Tuple>(tp) );
}

namespace tum {
  template < typename Tuple1, typename Tuple2, class = tum::enable_if_tuple_t<Tuple1,Tuple2> >
  constexpr auto inner_product ( Tuple1&& tp1, Tuple2&& tp2 ) {
    return std::apply( []( auto a, auto b ){return a+b;}, std::forward<Tuple1>(tp1) * std::forward<Tuple2>(tp2) );
  }

  template < typename Tuple, class = enable_if_tuple_t<Tuple> >
  constexpr auto abs_sq ( Tuple&& tp ) {
    return inner_product( std::forward<Tuple>(tp), std::forward<Tuple>(tp) );
  }

  template < typename Tuple, class = enable_if_tuple_t<Tuple> >
  constexpr auto abs ( Tuple&& tp ) {
    return std::sqrt( abs_sq( std::forward<Tuple>(tp) ) );
  }
}

#endif
