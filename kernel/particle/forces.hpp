#ifndef _PARTICLE_FORCES_HPP_
#define _PARTICLE_FORCES_HPP_

#include "apt/vec.hpp"
#include "particle/map.hpp"
#include "particle/array.hpp"
#include <vector>

namespace particle {
  // NOTE currently we assume all forces have at most one free param
  template < typename T, template < typename > class S >
  using force_t = void (*) ( typename array<T,S>::particle_type& ptc,
                             T dt,
                             const apt::Vec<T, S<T>::Dim>& E,
                             const apt::Vec<T, S<T>::Dim>& B,
                             T param0 );
}

namespace particle {
  template < typename T, template < typename > class S >
  struct Force {
    std::vector<force_t<T,S>> forces;
    std::vector<T> params;

    void add ( force_t<T,S> force, T param );

    void Register( species sp ) const;
    static void Unregister( species sp );
    static Force<T,S> Get( species sp ) noexcept;

  };
}

#define LORENTZ

#ifdef LORENTZ
#include <sstream>
#endif

namespace particle::force {
#ifdef LORENTZ
  extern std::ostringstream ostr;
#endif
  template < typename T, template < typename > class S, template < typename, template < typename > class > class Ptc_t >
  void lorentz( Ptc_t<T,S>& ptc, T dt, const apt::Vec<T, S<T>::Dim>& E, const apt::Vec<T, S<T>::Dim>& B, T q_times_w_gyro_unitB_over_m  );

  template < typename T, template < typename > class S, template < typename, template < typename > class > class Ptc_t >
  void lorentz_exact( Ptc_t<T,S>& ptc, T dt, const apt::Vec<T, S<T>::Dim>& E, const apt::Vec<T, S<T>::Dim>& B, T q_times_w_gyro_unitB_over_m  );
}


#endif
