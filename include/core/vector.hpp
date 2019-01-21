#ifndef  _VECTOR_HPP_
#define  _VECTOR_HPP_

#include <array>
#include <experimental/array> // for make_array
#include <cmath>
#include <tuple>
#include <type_traits>

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

template < typename T, std::size_t N >
struct Vec : public std::array<T,N> {
  template < typename U >
  constexpr Vec( std::array<U,N> x ) noexcept
    : std::array<T,N>( std::move(x) ) {}
  // TODO rule of five?

  // converting from reference array
  explicit constexpr Vec( const Vec<T&,N>& other ) noexcept {
    // TODO finish this
  }

};


template < typename T, std::size_t N >
struct Vec<T&, N> : public homogeneous_tuple<T&,N> {
  constexpr Vec( homogeneous_tuple<T&,N> x ) noexcept
    : homogeneous_tuple<T&,N>( std::move(x) ) {}

  // TODO double check the semantics of these
  Vec( const Vec& other ) = default;
  Vec( Vec&& other ) noexcept = default;
};

template < typename OStream, typename T, std::size_t N, std::size_t I=0 >
OStream& operator<< ( OStream& os, const Vec<T,N>& tp ) {
  if constexpr( I < N ) {
      os << std::get<I>(tp);
      if constexpr ( I < N-1 ) os << ", ";
      operator<< < OStream, T, N, I+1 > ( os, tp );
    }
  return os;
}

namespace std {
  // NOTE when a derived class is used as a template argument, there is just no conversion to its base class ever.
  template < std::size_t I, typename T, std::size_t N >
  struct tuple_element < I, Vec<T,N> > {
    using type = T;
  };

  template < typename T, std::size_t N >
  struct tuple_size < Vec<T,N> > {
    static constexpr auto value = N;
  };
}

namespace vec {
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

namespace vec {

  template < class Vec >
  using element_t = std::tuple_element_t<0, std::remove_reference_t<Vec> >;

  template < class Vec >
  inline constexpr auto size_v = std::tuple_size_v< std::remove_reference_t<Vec> >;

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

namespace vec :: per_dim {
  namespace impl {

    template < std::size_t D, typename Func, typename... Tuples >
    constexpr auto invoke_on_dim( const Func& f, Tuples&&... tpls ) noexcept {
      return f( std::get<D>(std::forward<Tuples>(tpls))... );
    }

    template < typename Func, typename... Tuples, std::size_t... D >
    constexpr auto make( const Func& f, std::index_sequence<D...>, Tuples&&... tpls ) noexcept {
      return vec::make( invoke_on_dim<D>( f, std::forward<Tuples>(tpls)... )... );
    }

    template < typename Func, typename... Tuples, std::size_t... D >
    constexpr auto tie( const Func& f, std::index_sequence<D...>, Tuples&... tpls ) noexcept {
      return vec::tie( invoke_on_dim<D>( f, tpls... )... );
    }

  }


  // f shouldn't depend on the order of each dim executed
  // TODO return Vec<T> or Vec<T&> ??? How to decide Vec<T>& vs Vec<T&>. Maybe just return Vec<T> or Vec<T&> ???
  // NOTE TODO currently per_dim_of requires Func to return something, and the return object of per_dim_of is always vec::made rather than vec::tied, so there is always a temporary object returned. Cannot use per_dim_of to generate reference
  template < std::size_t Ndims, typename Func, typename... Vectors >
  constexpr auto make( const Func& f, Vectors&&... vecs  ) noexcept {
    return impl::make( f, std::make_index_sequence<Ndims>{},
                       std::forward<Vectors>(vecs)... );
  }

  template < std::size_t Ndims, typename Func, typename... Vectors >
  constexpr auto tie( const Func& f, Vectors&... vecs  ) noexcept {
    return impl::tie( f, std::make_index_sequence<Ndims>{}, vecs... );
  }

}

namespace vec {
  // NOTE foreach has no return type as opposed to those in vec::per_dim
  template < std::size_t Begin, std::size_t End, typename Func, typename... Vectors >
  constexpr void foreach( const Func& f, Vectors&&... vecs  ) noexcept {
    static_assert( Begin <= End );
    if constexpr ( Begin == End ) return;
    else {
      f( std::get<Begin>(std::forward<Vectors>(vecs))... );
      return foreach<Begin+1, End>( f, std::forward<Vectors>(vecs) );
    }
  }
}


namespace vec::functional {
  template < typename Vector, typename Vec_or_Num, typename BinaryOp >
  constexpr auto adapt_binary_by_make( const BinaryOp& op ) noexcept {
    static_assert( is_vec_v<Vector>, "First argument must be of Vec type." );
    if constexpr ( is_vec_v<Vec_or_Num> ) {
      constexpr auto N = std::min( size_v<Vector>, size_v<Vec_or_Num> );
      return [&op]( const Vector& v1, const Vec_or_Num& v2 ) noexcept {
               return vec::per_dim::make<N> ( op, v1, v2 );};
    } else  {
      static_assert( std::is_arithmetic_v<Vec_or_Num> );
      using T = element_t<Vector>;
      constexpr auto N = size_v<Vector>;
      return [&op]( const Vector& v, const Vec_or_Num& s ) {
               auto&& f = [ &op, &rhs = s ] ( const T& lhs ) noexcept {
                            return op(lhs, rhs); };
               return vec::per_dim::make<N> ( std::move(f), v );
             };
    }
  }

  template < typename Vector, typename Vec_or_Num, typename BinaryOp >
  constexpr auto adapt_binary_by_tie( const BinaryOp& op ) noexcept {
    static_assert( is_vec_v<Vector>, "First argument must be of Vec type." );
    if constexpr ( is_vec_v<Vec_or_Num> ) {
      constexpr auto N = std::min( size_v<Vector>, size_v<Vec_or_Num> );
      return [&op]( Vector& v1, const Vec_or_Num& v2 ) noexcept {
               return vec::per_dim::tie<N> ( op, v1, v2 );};
    } else  {
      static_assert( std::is_arithmetic_v<Vec_or_Num> );
      using T = element_t<Vector>;
      constexpr auto N = size_v<Vector>;
      return [&op]( Vector& v, const Vec_or_Num& s ) noexcept {
               auto&& f = [ &op, &rhs = s ] ( T& lhs ) noexcept {
                            return op(lhs, rhs); };
               return vec::per_dim::tie<N> ( std::move(f), v );
             };
    }
  }

  template < class T, typename MemGet_F >
  constexpr auto mem_get( const MemGet_F& get_mem ) {
    return [&get_mem] ( T& obj ) noexcept {
             return vec::per_dim::tie<size_v<T>>( get_mem, obj );
           };
  }

}

namespace {
// NOTE somehow generic lambda works here!
#define vec_def_binary_op( _OP_, _MODE_ )                                   \
  template < typename T, std::size_t N, typename U >                    \
  constexpr auto operator _OP_ ( const Vec<T,N>& v1, const U& v2 ) {    \
    return vec::functional::adapt_binary_##_MODE_<Vec<T,N>, U>          \
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
    return vec::per_dim::tie<vec::size_v<_CLASS_>>    \
      ( []( auto&& elm ){ return elm._MEM_;}, obj );  \
  }
}


namespace vec::numeric {
  template < typename T1, typename T2, std::size_t N >
  constexpr auto dot ( const Vec<T1,N>& v1, const Vec<T2,N>& v2 ) {
    auto&& tmp = vec::per_dim::make<N>( []( auto a, auto b ) {return a*b;}, std::move(v1), std::move(v2) );
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
    return vec::make( get<1>(v1) * get<2>(v2) - get<2>(v1) * get<1>(v2),
                      get<2>(v1) * get<0>(v2) - get<0>(v1) * get<2>(v2),
                      get<0>(v1) * get<1>(v2) - get<1>(v1) * get<0>(v2) );
  }
}

namespace vec {
  namespace impl {
    template < typename From, typename To >
    struct copy_ref {
      using type = std::conditional_t< std::is_reference_v<From>, To&, To >;
    };

    template < typename From, typename To >
    struct copy_const {
      using type = std::conditional_t< std::is_const_v<From>, const To, To >;
    };

  }

  template < typename From, typename To >
  using copy_ref_t = typename copy_ref<From,To>::type;

  template < typename From, typename To >
  using copy_const_t = typename copy_const<From,To>::type;

  template < typename From, typename To >
  using copy_constref_t = typename copy_const_t<From, copy_ref_t<From, To>>;

  // NOTE: somehow I don't want to use std::decay. This template will be provided in std in C++20
  template < typename T >
  using remove_cvref_t = std::remove_cv_t<std::remove_reference_t<T>>;
}

#endif
