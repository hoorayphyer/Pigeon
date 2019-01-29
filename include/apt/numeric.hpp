#ifndef  _APT_NUMERIC_HPP_
#define  _APT_NUMERIC_HPP_

namespace apt {
  template < typename T1, typename T2, std::size_t N >
  constexpr auto dot ( const Vec<T1,N>& v1, const Vec<T2,N>& v2 ) {
    auto&& tmp = apt::per_dim::make<N>( []( auto a, auto b ) {return a*b;}, std::move(v1), std::move(v2) );
    return std::apply( [](auto&&... args){
                         return (... + std::forward<decltype(args)>(args) );
                       }, std::move(tmp) );
  }

  template < typename T, std::size_t N >
  constexpr auto abs_sq ( const Vec<T,N>& v ) {
    static_assert( N < 4, "not implemented" );
    if constexpr ( N == 1 ) return std::abs( std::get<0>(v) );
    else return std::apply( std::hypot, v );
  }

  template < typename T, std::size_t N >
  constexpr auto abs ( const Vec<T,N>& v ) {
    return std::sqrt( abs_sq(v) );
  }

  template < typename T1, typename T2, std::size_t N >
  constexpr auto cross ( const Vec<T1,N>& v1, const Vec<T2,N>& v2 ) {
    static_assert( N == 3, "not implemented" );
    using std::get;
    return apt::make( get<1>(v1) * get<2>(v2) - get<2>(v1) * get<1>(v2),
                      get<2>(v1) * get<0>(v2) - get<0>(v1) * get<2>(v2),
                      get<0>(v1) * get<1>(v2) - get<1>(v1) * get<0>(v2) );
  }
}

#endif
