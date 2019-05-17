#ifndef _MISC_HPP_
#define _MISC_HPP_
#include <cmath>

namespace misc {
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
    const auto origin = mesh.origin();
    for ( int i = 0; i < DField; ++i ) {
      for ( auto I : apt::Block( mesh.extent() ) ) {
        I += origin;
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
}

#endif
