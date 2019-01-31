#ifndef  _APT_NUMERIC_HPP_
#define  _APT_NUMERIC_HPP_
#include "apt/vec.hpp"

// TODO use expression template to optimize
namespace apt::functional {
  template < typename Vector, typename Vec_or_Num, typename BinaryOp >
  constexpr auto adapt_binary_by_make( const BinaryOp& op ) noexcept {
    static_assert( is_vec_v<Vector>, "First argument must be of Vec type." );
    if constexpr ( is_vec_v<Vec_or_Num> ) {
      constexpr auto N = std::min( size_v<Vector>, size_v<Vec_or_Num> );
      return [&op]( const Vector& v1, const Vec_or_Num& v2 ) noexcept {
               return apt::per_dim::make<N> ( op, v1, v2 );};
    } else  {
      static_assert( std::is_arithmetic_v<Vec_or_Num> );
      using T = element_t<Vector>;
      constexpr auto N = size_v<Vector>;
      return [&op]( const Vector& v, const Vec_or_Num& s ) {
               auto&& f = [ &op, &rhs = s ] ( const T& lhs ) noexcept {
                            return op(lhs, rhs); };
               return apt::per_dim::make<N> ( std::move(f), v );
             };
    }
  }

  template < typename Vector, typename Vec_or_Num, typename BinaryOp >
  constexpr auto adapt_binary_by_tie( const BinaryOp& op ) noexcept {
    static_assert( is_vec_v<Vector>, "First argument must be of Vec type." );
    if constexpr ( is_vec_v<Vec_or_Num> ) {
      constexpr auto N = std::min( size_v<Vector>, size_v<Vec_or_Num> );
      return [&op]( Vector& v1, const Vec_or_Num& v2 ) noexcept {
               return apt::per_dim::tie<N> ( op, v1, v2 );};
    } else  {
      static_assert( std::is_arithmetic_v<Vec_or_Num> );
      using T = element_t<Vector>;
      constexpr auto N = size_v<Vector>;
      return [&op]( Vector& v, const Vec_or_Num& s ) noexcept {
               auto&& f = [ &op, &rhs = s ] ( T& lhs ) noexcept {
                            return op(lhs, rhs); };
               return apt::per_dim::tie<N> ( std::move(f), v );
             };
    }
  }

  template < class T, typename MemGet_F >
  constexpr auto mem_get( const MemGet_F& get_mem ) {
    return [&get_mem] ( T& obj ) noexcept {
             return apt::per_dim::tie<size_v<T>>( get_mem, obj );
           };
  }

}

namespace {
// NOTE somehow generic lambda works here!
#define vec_def_binary_op( _OP_, _MODE_ )                                   \
  template < typename T, std::size_t N, typename U >                    \
  constexpr auto operator _OP_ ( const Vec<T,N>& v1, const U& v2 ) {    \
    return apt::functional::adapt_binary_##_MODE_<Vec<T,N>, U>          \
      ([]( auto&& a, auto&& b ) { return a _OP_ b; }) ( v1, v2 );       \
  }                                                                     \

  vec_def_binary_op( + , by_make );
  vec_def_binary_op( - , by_make );
  vec_def_binary_op( * , by_make );
  vec_def_binary_op( / , by_make );

  vec_def_binary_op( += , by_tie );
  vec_def_binary_op( -= , by_tie );
  vec_def_binary_op( *= , by_tie );
  vec_def_binary_op( /= , by_tie );

  template < typename T, std::size_t N, typename U,
             class = std::enable_if_t< std::is_arithmetic_v<U>, int > >
  constexpr auto operator* ( U val, const Vec<T,N>& v ) {
    return v * val;
  }
#undef vec_def_binary_op

#define vec_def_member_getter( _CLASS_, _MEM_)        \
  constexpr auto _MEM_( _CLASS_& obj ) noexcept {     \
    return apt::per_dim::tie<apt::size_v<_CLASS_>>    \
      ( []( auto&& elm ){ return elm._MEM_;}, obj );  \
  }
}

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
