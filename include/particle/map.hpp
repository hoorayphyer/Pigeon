#ifndef  _PARTICLE_MAP_HPP_
#define  _PARTICLE_MAP_HPP_

#include "particle/species_predef.hpp"
#include <unordered_map>
#include <initializer_list>

namespace particle {

  template < typename Val >
  struct map {
  private:
    std::unordered_map<species, Val> _data;
  public:
    map () = default;

    map ( std::initializer_list<species> list ) {
      for ( auto sp : list ) _data.emplace( sp, {} );
    }

    template < typename U >
    inline void add_if_not( species sp ) noexcept {
      if ( _data.find(sp) != _data.end() )
        _data.emplace( sp, {} );
    }

    inline Val& operator[] ( species sp ) noexcept {
      return _data.at(sp);
    }

    inline const Val& operator[] ( species sp ) const noexcept {
      return _data.at(sp);
    }

    inline bool has( species sp ) const noexcept {
      return _data.find(sp) != _data.end();
    }

    // TODO
    auto begin() noexcept {
      
    }

    auto end() noexcept {
      
    }
  };
}

#endif
