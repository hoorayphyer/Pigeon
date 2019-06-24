#ifndef _PARTICLE_SORTER_HPP_
#define _PARTICLE_SORTER_HPP_
#include "particle/array.hpp"

namespace particle {
  // TODO the function now only erases empty particles. Also add sorting particles by position
  template < typename T, template < typename > class Spec >
  void sort(  array<T,Spec>& ptcs ) {
    int head = 0;
    int tail = ptcs.size() - 1;
    while( head < tail ) {
      if ( ptcs[head].is(flag::exist) ) {
        ++head;
      } else {
        if ( !ptcs[tail].is(flag::exist) ) {
          --tail;
        } else {
          ptcs[head] = std::move(ptcs[tail]);
          ++head;
          --tail;
        }
      }
    }
    // when exit head = new size
    ptcs.resize(head);
  }
}

#endif
