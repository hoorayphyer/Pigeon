#ifndef _APT_PRINT_VEC_HPP_
#define _APT_PRINT_VEC_HPP_

#include "apt/vec_expression.hpp"
#include <string>

namespace std {
  template < typename E, typename T, bool L >
  string to_string( const apt::VecExpression<E,T,L>& vec ) {
    constexpr int N = E::size;
    string res = "( " + to_string( std::get<0>(vec) );
    apt::foreach<1,N>
      ( [&]( auto x ) {
          res += ", " + to_string(x);
        }, vec );
    res += " )";
    return res;
  }
}

template < typename OStream, typename E, typename T, bool L >
OStream& operator<< ( OStream& os, const apt::VecExpression<E,T,L>& vec ) {
  constexpr int N = E::size;
  os << "( " << std::get<0>(vec);
  apt::foreach<1,N>
    ( [&]( auto x ) {
        os << ", " << x;
      }, vec );
  os << " )";
  return os;
}

#endif
