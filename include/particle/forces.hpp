#ifndef _PARTICLE_FORCES_HPP_
#define _PARTICLE_FORCES_HPP_

#include "apt/vec.hpp"
#include "particle/map.hpp"
#include "particle/virtual_particle.hpp"
#include <vector>

namespace particle::force {
  template < typename T, template < typename > class PtcSpecs >
  using Vec = apt::Vec<T, PtcSpecs<T>::Dim>;

  // NOTE currently we assume all forces have at most one free param
  template < typename T, template < typename > class PtcSpecs, template < typename, template < typename > class > class Ptc_t >
  using force_t = void (*) ( Ptc_t<T,PtcSpecs>& ptc,
                             T dt,
                             const Vec<T,PtcSpecs>& E,
                             const Vec<T,PtcSpecs>& B,
                             T param0 );
}

namespace particle {
  template < typename T, template < typename > class PtcSpecs, template < typename, template < typename > class > class Ptc_t >
  struct ForceGen;

  namespace force {
    template < typename T, template < typename > class PtcSpecs, template < typename, template < typename > class > class Ptc_t >
    struct Force {
    private:

      std::vector<force_t<T,PtcSpecs,Ptc_t>> _forces;
      std::vector<T> _params;

    public:
      friend class ForceGen<T,PtcSpecs, Ptc_t>;
      // TODO turn these into unique_ptr and use deepcopy
      void add ( force::force_t<T,PtcSpecs,Ptc_t> force, T param ) {
        _forces.push_back(force);
        _params.push_back(param);
      }
    };
  }


  template < typename T, template < typename > class PtcSpecs, template < typename, template < typename > class > class Ptc_t >
  struct ForceGen {
  private:
    map<force::Force<T,PtcSpecs,Ptc_t>> _force_map;

  public:
    void Register( species sp, force::Force<T,PtcSpecs,Ptc_t> force ) {
      _force_map[sp] = std::move(force);
    }

    void Unregister( species sp ) {
      _force_map.erase(sp);
    }

    inline auto operator() ( species sp ) noexcept {
      // if sp doesn't exist, newly created item in the maps will automatically have zero elements
      return [&forces = _force_map[sp]._forces,
              &params = _force_map[sp]._params]
        ( auto& ptc, auto&&... args ) {
          for ( int i = 0; i < forces.size(); ++i ) {
            (forces[i])( ptc, std::forward<decltype(args)>(args)..., params[i] );
          }
        };
    }

  };
}

namespace particle::force {
  template < typename T, template < typename > class PtcSpecs, template < typename, template < typename > class > class Ptc_t >
  void lorentz( Ptc_t<T,PtcSpecs>& ptc, T dt, const Vec<T,PtcSpecs>& E, const Vec<T,PtcSpecs>& B, T q_over_m  );
}


#endif
