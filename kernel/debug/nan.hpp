#ifndef _DEBUGGER_NAN_HPP_
#define _DEBUGGER_NAN_HPP_

#include "apt/block.hpp"
#include "apt/range.hpp"
#include <cmath>

namespace debug {
  template < typename T, int DField, int DGrid, template < typename, int, int > class Field >
  inline int num_nan( const Field<T,DField,DGrid>& f ) {
    int res = 0;
    for ( int i = 0; i < DField; ++i ) {
      for( const auto& x : f[i].data() )
        res += std::isnan(x);
    }
    return res;
  }

  template < typename T, int DField, int DGrid, template < typename, int, int > class Field, class Logger >
  inline int show_nan( const Field<T,DField,DGrid>& f, Logger& log ) {
    int num = 0;
    const auto& mesh = f.mesh();
    for ( int i = 0; i < DField; ++i ) {
      for ( auto I : apt::Block<DGrid>( apt::range::far_begin(mesh.range()), apt::range::far_end(mesh.range()) ) ) {
        if ( std::isnan(f[i](I)) ) {
          log % "  Found NAN! comp=" << i << ", indices=" << I << std::endl;
          ++num;
          if ( num >= 20 ) {
            log % "  More than " << num << " NANs found at comp=" << i << ", break" << std::endl;
            break;
          }
        }
      }
    }
    return num;
  }

  template < typename T, template < typename > class Specs, template < typename, template < typename > class > class Ptc >
  inline int has_nan( const Ptc<T,Specs>& ptc ) {
    int res = 0;
    constexpr auto DPtc = Specs<T>::Dim;
    for ( int i = 0; i < DPtc; ++i ) {
      if ( std::isnan(ptc.q()[i]) ) res |= ( 1u << i);
      if ( std::isnan(ptc.p()[i]) ) res |= ( 1u << (i+DPtc) );
    }
    if ( std::isnan(ptc.state() )) res |= ( 1u << (2 * DPtc) );
    return res;
  }
}

#endif
