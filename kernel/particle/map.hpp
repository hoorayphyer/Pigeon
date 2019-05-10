#ifndef  _PARTICLE_MAP_HPP_
#define  _PARTICLE_MAP_HPP_

#include "particle/species_predef.hpp"
#include <unordered_map>

namespace particle {

  template < typename Val >
  struct map : std::unordered_map<particle::species, Val> {
    using std::unordered_map<particle::species, Val>::unordered_map;

    inline bool has( species sp ) const noexcept {
      return this->find(sp) != this->end();
    }
  };
}

#endif
