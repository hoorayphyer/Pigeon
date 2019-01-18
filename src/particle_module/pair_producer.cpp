#include "particle_module/pair_producer.hpp"

namespace particle {

  inline Real SamplePhotonEnergy(const Real &gamma, const Real &Rc) {
    return _pane.E_ph;
  }

  inline bool IsProduceCurvPhoton(Real dt, const Real &gamma, const float &Rc, Real randnum ) {
    Real Gamma_thr = _pane.K_thr *  std::cbrt(Rc);
    return ( gamma > _pane.gamma_off ) && ( gamma > Gamma_thr )
      && ( randnum  < _pane.K_curv_em_rate * gamma * dt / Rc );
  }

  template <typename Ptc >
  void instant_produce_pairs( Real dt, Rng& rng, Ptc& ptc, Ptc& electron, Ptc& positron ) {
    namespace vn = vec::numeric;

    if( _pane.is_ignore_pair_production( ptc.q ) ) continue;

    Real pimag = std::sqrt( vn::abs_sq( ptc.p ) );
    Real gamma_i = std::sqrt( 1.0 + pimag * pimag );

    // check if curvature photons capable of pair producing can be created.
    if ( !IsProduceCurvPhoton(dt, gamma_i, Rc, rng.uniform()) ) continue;

    Real E_ph = SamplePhotonEnergy( gamma_i, Rc );
    Real gamma_f = gamma_i - E_ph;
    Real pfmag = std::sqrt( gamma_f * gamma_f - 1.0 );
    Real del_pmag = pimag - pfmag;// note del_pmag > 0

    // create pairs
    Real InstantPairCreateImpl( const Real &del_gamma, const Real &del_p ) {
      // equally partition the energy, allowing momentum to be not conserved.
      Real gamma_sec = del_gamma / 2.0;
      return std::sqrt( gamma_sec * gamma_sec - 1 );
      pp_ptron = _etron = p_ptron;
    }
    Real p_ptron = InstantPairCreateImpl(E_ph, del_pmag);;
    Real p_etron = p_ptron;

    // append electron and positron
    electron.q = ptc.q;
    electron.p = ptc.p * (p_etron / pimag);
    electron.set<flag::secondary>();

    positron.q = electron.q;
    positron.p = electron.p;
    positron.set<flag::secondary>();

    // primary particle loses energy to gamma rays
    p *= ( pfmag / pimag );

    // // track the pair creation event if the primary is tracked
    // if ( ParticleTracker::IsTracked(ptc.track_id) ) {
    //   e_tracker.Track( ptc_e );
    //   p_tracker.Track( ptc_p );
    //   auto parentType = particles.Attributes().ptcType;
    //   pc_tracker.LogThisEvent( timestep, parentType, ptc.track_id, ParticleType::ELECTRON, ptc_e.track_id, ParticleType::POSITRON, ptc_p.track_id );
    // }

    // TODO
    // electrons.Append( std::move(ptc_e) );
    // positrons.Append( std::move(ptc_p) );


    // TODO register this pair creation event
    // data.pairCreationEvents.data() [cell] += 1.0;

  }

  inline float f_x ( Real x ) {
    // distribution of x*exp(-x^2/2), which peaks at x = 1.
    return std::sqrt( -2.0 * std::log(x) );
  }

  inline bool IsMagneticConvert( Real dt, const Vec3<Real>& Bvec, Real randnum ) {
    return Bvec.dot(Bvec) > _pane.B_magconv * _pane.B_magconv ? ( randnum < dt / _pane.l_magconv ) : false;
  }

  template <typename Ptc >
  void photon_produce_pairs( Real dt, Ptc& photon, Ptc& electron, Ptc& positron, Rng& rng ) {
    Vec3<Real> BVec ( data.Bfield.data(0)[cell], data.Bfield.data(1)[cell], data.Bfield.data(2)[cell]);
    if ( !IsMagneticConvert( dt, BVec, rng.uniform() ) && photons.PtcData()[idx].path_left > 0.0 ) continue;

    // ratio = pmag_ptc / E_ph. Here the pair production process is assumed to conserve energy but not momentum.
    Real&& ratio = std::sqrt( 0.25 - 1.0 / vn::abs_sq(mem::p(photon)) );

    // append electron and positron

    electron.q = photon.q;
    electron.p = photon.p * std::move(ratio);
    electron.set<flag::secondary>();

    positron.q = electron.q;
    positron.p = electron.p;
    positron.set<flag::secondary>();

    // // only track the event if the parent is tracked
    // if ( ParticleTracker::IsTracked(ptc.track_id) ) {
    //   auto& last_e = electrons.PtcData()[ electrons.Number() - 1 ];
    //   auto& last_p = positrons.PtcData()[ positrons.Number() - 1 ];
    //   e_tracker.Track( last_e );
    //   p_tracker.Track( last_p );
    //   pc_tracker.LogThisEvent( timestep, ParticleType::PHOTON, ptc.track_id, ParticleType::ELECTRON, last_e.track_id, ParticleType::POSITRON, last_p.track_id );
    // }

    // // register this pair creation event
    // data.pairCreationEvents.data() [cell] += 1.0;

    // erase this photon
    photon.set<flag::nonexist>();
  }

  template <typename Ptc >
  void produce_photons( Ptc& ptc, Ptc& photon, Rng& rng ) {
    // carry out production
    const auto& Rc = ptc.Rc;

    if( _pane.is_ignore_pair_production( ptc.q ) ) continue;

    Real pimag = vn::abs( ptc.p);
    Real gamma_i = std::sqrt( 1.0 + pimag * pimag );

    // check if curvature photons capable of pair producing can be created.
    if ( !IsProduceCurvPhoton(dt, gamma_i, Rc, rng.uniform()) ) continue;

    Real E_ph = SamplePhotonEnergy( gamma_i, Rc );
    Real gamma_f = std::max(1.0, gamma_i - E_ph);
    Real pfmag = std::sqrt( gamma_f * gamma_f - 1.0 );

    // create photons
    float path_left = _pane.l_coll * f_x( rng.uniform() );
    photon.q = ptc.q;
    photon.p = ptc.p * (E_ph / pimag);
    ph_new.path_left = path_left;

    // primary particle loses energy to gamma rays
    ptc.p *= ( pfmag / pimag );

    // only track the event if the parent is tracked
    // if ( ParticleTracker::IsTracked(ptc.track_id) ) {
    //   auto& last_ph = photons.PtcData()[ photons.Number() - 1 ];
    //   ph_tracker.Track( last_ph );
    //   auto parentType = particles.Attributes().ptcType;
    //   pc_tracker.LogThisEvent( timestep, parentType, ptc.track_id, ParticleType::PHOTON, last_ph.track_id );
    // }

  }

}
