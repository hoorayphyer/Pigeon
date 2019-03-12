#ifndef _APT_PRINT_HPP_
#define _APT_PRINT_HPP_

#include "apt/array.hpp"
#include "apt/vec_expression.hpp"
#include <string>

namespace std {
  template < typename OStream, typename T, int N >
  string to_string( const apt::array<T,N>& c ) {
    string res = "( " + to_string( c[0] );
    for ( int i = 1; i < N; ++i )
      res += (", " + to_string(c[i]));
    res += " )";
    return res;
  }

  template < typename OStream, typename E >
  string to_string( const apt::VecExpression<E>& c ) {
    string res = "( " + to_string( c[0] );
    for ( int i = 1; i < E::NDim; ++i )
      res += (", " + to_string(c[i]));
    res += " )";
    return res;
  }
}

template < typename OStream, typename T, int N >
OStream& operator<< ( OStream& os, const apt::array<T,N>& c ) {
  os << "( " << c[0];
  for ( int i = 1; i < N; ++i )
    os << ", " << c[i];
  os << " )";
  return os;
}

template < typename OStream, typename E >
OStream& operator<< ( OStream& os, const apt::VecExpression<E>& c ) {
  os << "( " << c[0];
  for ( int i = 1; i < E::NDim; ++i )
    os << ", " << c[i];
  os << " )";
  return os;
}

#endif
