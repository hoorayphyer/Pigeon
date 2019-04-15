#ifndef _PARTICLE_FORCES_HPP_
#define _PARTICLE_FORCES_HPP_

#include "apt/vec.hpp"
#include "particle/map.hpp"
#include <vector>

namespace particle::force {
  // TypeStruct for forces
  namespace ts {
    template < class Ptc >
    using Real = typename Ptc::vec_type::element_type;

    template < class Ptc >
    using Vec = apt::Vec<Real<Ptc>,Ptc::NDim>;
  };

  // NOTE currently we assume all forces have at most one free param
  template < class Ptc >
  using force_t = void (*) ( Ptc& ptc,
                             ts::Real<Ptc> dt,
                             const ts::Vec<Ptc>& E,
                             const ts::Vec<Ptc>& B,
                             ts::Real<Ptc> param0 );
}

namespace particle {
  template < class Ptc >
  struct ForceGen;

  namespace force {
    template < class Ptc >
    struct Force {
    private:
      using Real = typename Ptc::vec_type::element_type;

      std::vector<force::force_t<Ptc>> _forces;
      std::vector<Real> _params;

    public:
      friend class ForceGen<Ptc>;
      // TODO turn these into unique_ptr and use deepcopy
      void add ( force::force_t<Ptc> force, Real param ) {
        _forces.push_back(force);
        _params.push_back(param);
      }
    };
  }


  template < class Ptc >
  struct ForceGen {
  private:
    using Real = typename Ptc::vec_type::element_type;
    map<force::Force<Ptc>> _force_map;

  public:
    void Register( species sp, force::Force<Ptc> force ) {
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
  template < typename Ptc >
  void lorentz( Ptc& ptc, ts::Real<Ptc> dt, const ts::Vec<Ptc>& E, const ts::Vec<Ptc>& B, ts::Real<Ptc> q_over_m  );

  // when B is strong enough, damp the perpendicular component of momentum
  template < typename Ptc >
  void landau0( Ptc& ptc, ts::Real<Ptc> dt, const ts::Vec<Ptc>& E, const ts::Vec<Ptc>& B, ts::Real<Ptc> B2_thr  );
}


#endif
