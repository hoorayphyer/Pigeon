#ifndef _PARTICLE_FORCES_HPP_
#define _PARTICLE_FORCES_HPP_

#include <string>
#include "apt/vec.hpp"

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

namespace particle::force {
  using id_t = std::string;

  template < class Ptc >
  struct specs {
    id_t id;
    ts::Real<Ptc> param0;
  };

  template < class Ptc >
  struct Factory {
    static void Register( const id_t& id, force_t<Ptc> force );
    static void Unregister( const id_t& id );
    static force_t<Ptc> create( const id_t& id );
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
