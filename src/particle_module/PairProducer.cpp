#include "PairProducer.h"
#include "AperData.h"
#include "Logger.h"
#include "InfoCollector.h"
#include "Dashboard.h"
#include "Rng.h"

PairProducer_Instant::PairProducer_Instant( const DBPane_PairProduction& pane ) : _pane(pane) {}

inline Scalar PairProducer_Instant::SamplePhotonEnergy(const Scalar &gamma, const float &Rc) const {
  return _pane.E_ph;
}

inline bool PairProducer_Instant::IsProduceCurvPhoton(Scalar dt, const Scalar &gamma, const float &Rc, Scalar randnum ) const {
  Scalar Gamma_thr = _pane.K_thr *  std::cbrt(Rc);
  return ( gamma > _pane.gamma_off ) && ( gamma > Gamma_thr )
    && ( randnum  < _pane.K_curv_em_rate * gamma * dt / Rc );
}

void InstantPairCreateImpl(Scalar &p_ptron, Scalar &p_etron, const Scalar &del_gamma, const Scalar &del_p) {
  // equally partition the energy, allowing momentum to be not conserved.
  Scalar gamma_sec = del_gamma / 2.0;
  p_ptron = std::sqrt( gamma_sec * gamma_sec - 1 );
  p_etron = p_ptron;
}

void PairProducer_Instant::ProducePairs( int timestep, AperData& data, Scalar dt, const Grid &grid, Rng& rng ) const {
  int num = 0;
  auto& electrons = data.particles.at(ParticleType::ELECTRON);
  auto& positrons = data.particles.at(ParticleType::POSITRON);
  auto& e_tracker = electrons.Tracker();
  auto& p_tracker = positrons.Tracker();
  auto& pc_tracker = data.pairCreationTracker;
  int pairs_flag = 0;
  set_bit( pairs_flag, ParticleFlag::secondary );
  // carry out production

  for ( auto& elm : data.particles ) {
    auto& particles = elm.second;
    if (!particles.Attributes().isRadiate) continue;
    for ( unsigned int idx = 0; idx < particles.Number(); ++idx ) {
      if( particles.IsEmpty(idx) ) continue;

      QCarrier& ptc = particles.PtcData()[idx];
      const float& Rc = ptc.Rc;
      Vec3<MOM_TYPE>& p = ptc.p;
      const int& cell = ptc.cell;
      const Vec3<POS_TYPE>& x = ptc.x;

      Vec3<Scalar> pos = grid.pos_particle( cell, x );
      if( _pane.is_ignore_pair_production( pos.x, pos.y, pos.z ) ) continue;

      MOM_TYPE pimag = std::sqrt( p.dot(p) );
      Scalar gamma_i = std::sqrt( 1.0 + pimag * pimag );

      // check if curvature photons capable of pair producing can be created.
      if ( !IsProduceCurvPhoton(dt, gamma_i, Rc, rng.uniform()) ) continue;

      Scalar E_ph = SamplePhotonEnergy( gamma_i, Rc );
      Scalar gamma_f = gamma_i - E_ph;
      MOM_TYPE pfmag = std::sqrt( gamma_f * gamma_f - 1.0 );
      MOM_TYPE del_pmag = pimag - pfmag;// note del_pmag > 0

      // create pairs
      Scalar p_ptron = 0.0;
      Scalar p_etron = 0.0;
      InstantPairCreateImpl(p_ptron, p_etron, E_ph, del_pmag);
      // append electron and positron
      QCarrier ptc_e;
      ptc_e.x = x, ptc_e.p = p * (p_etron / pimag), ptc_e.cell = cell, ptc_e.flag = pairs_flag;

      QCarrier ptc_p;
      ptc_p.x = x, ptc_p.p = p * (p_ptron / pimag), ptc_p.cell = cell, ptc_p.flag = pairs_flag;

      // track the pair creation event if the primary is tracked
      if ( ParticleTracker::IsTracked(ptc.track_id) ) {
        e_tracker.Track( ptc_e );
        p_tracker.Track( ptc_p );
        auto parentType = particles.Attributes().ptcType;
        pc_tracker.LogThisEvent( timestep, parentType, ptc.track_id, ParticleType::ELECTRON, ptc_e.track_id, ParticleType::POSITRON, ptc_p.track_id );
      }

      electrons.Append( std::move(ptc_e) );
      positrons.Append( std::move(ptc_p) );

      num++;

      // register this pair creation event
      data.pairCreationEvents.data() [cell] += 1.0;

      // primary particle loses energy to gamma rays
      p *= ( pfmag / pimag );
    }
  }


  Logger::print( 0, num, "pairs produced in this time step!" );
  InfoCollector::Instance().contents.num_createdPairs = num;

}

PairProducer_Photon::PairProducer_Photon( const DBPane_PairProduction& pane ) : PairProducer_Instant( pane ) {}

inline float f_x ( Scalar x ) {
  // distribution of x*exp(-x^2/2), which peaks at x = 1.
  return std::sqrt( -2.0 * std::log(x) );
}

inline bool PairProducer_Photon::IsMagneticConvert( Scalar dt, const Vec3<Scalar>& Bvec, Scalar randnum ) const {
  return Bvec.dot(Bvec) > _pane.B_magconv * _pane.B_magconv ? ( randnum < dt / _pane.l_magconv ) : false;
}

void PairProducer_Photon::ProducePairs( int timestep, AperData& data, Scalar dt, const Grid &grid, Rng& rng) const {
  int numPairs = 0;
  auto& photons = data.photons;
  auto& electrons = data.particles.at(ParticleType::ELECTRON);
  auto& positrons = data.particles.at(ParticleType::POSITRON);
  auto& e_tracker = electrons.Tracker();
  auto& p_tracker = positrons.Tracker();
  auto& pc_tracker = data.pairCreationTracker;
  int pairs_flag = 0;
  set_bit( pairs_flag, ParticleFlag::secondary );
  for ( unsigned int idx = 0; idx < photons.Number(); ++idx ) {
    if ( photons.IsEmpty( idx ) ) continue;
    const QNeutral& ptc = photons.PtcData()[idx];
    const int& cell = ptc.cell;
    Vec3<Scalar> BVec ( data.Bfield.data(0)[cell], data.Bfield.data(1)[cell], data.Bfield.data(2)[cell]);
    if ( !IsMagneticConvert( dt, BVec, rng.uniform() ) && photons.PtcData()[idx].path_left > 0.0 ) continue;

    const Vec3<POS_TYPE>& x = ptc.x;
    const auto& E_ph = ptc.E;
    const Vec3<PH_MOM_TYPE>& p_ph = ptc.p;

    // ratio = pmag_ptc / E_ph. Here the pair production process is assumed to conserve energy but not momentum.
    Scalar ratio = std::sqrt( 0.25 - 1.0 / ( E_ph * E_ph ) );
    Vec3<MOM_TYPE> p_ptc( p_ph * ratio );

    // append electron and positron
    QCarrier ptc_2nd;
    ptc_2nd.x = x, ptc_2nd.p = p_ptc, ptc_2nd.cell = cell, ptc_2nd.flag = pairs_flag;
    electrons.Append( ptc_2nd );
    positrons.Append( ptc_2nd );
    // only track the event if the parent is tracked
    if ( ParticleTracker::IsTracked(ptc.track_id) ) {
      auto& last_e = electrons.PtcData()[ electrons.Number() - 1 ];
      auto& last_p = positrons.PtcData()[ positrons.Number() - 1 ];
      e_tracker.Track( last_e );
      p_tracker.Track( last_p );
      pc_tracker.LogThisEvent( timestep, ParticleType::PHOTON, ptc.track_id, ParticleType::ELECTRON, last_e.track_id, ParticleType::POSITRON, last_p.track_id );
    }

    // register this pair creation event
    data.pairCreationEvents.data() [cell] += 1.0;

    // erase this photon
    photons.Erase( idx, 1 );

    numPairs++;
  }
  Logger::print( 0, numPairs, "pairs produced in this time step!" );
  InfoCollector::Instance().contents.num_createdPairs = numPairs;
}

void PairProducer_Photon::ProducePhotons( int timestep, AperData& data, Scalar dt, const Grid &grid, Rng& rng ) const {
  int num = 0;
  auto& photons = data.photons;
  auto& ph_tracker = photons.Tracker();
  auto& pc_tracker = data.pairCreationTracker;
  int photon_flag = 0;
  // carry out production
  for ( auto& elm : data.particles ) {
    auto& particles = elm.second;
    if (!particles.Attributes().isRadiate) continue;
    for ( unsigned int idx = 0; idx < particles.Number(); ++idx ) {
      if( particles.IsEmpty(idx) ) continue;

      QCarrier& ptc = particles.PtcData()[idx];
      const auto& Rc = ptc.Rc;
      Vec3<MOM_TYPE>& p = ptc.p;
      const int& cell = ptc.cell;
      const Vec3<POS_TYPE>& x = ptc.x;

      Vec3<Scalar> pos = grid.pos_particle( cell, x );
      if( _pane.is_ignore_pair_production( pos.x, pos.y, pos.z ) ) continue;

      MOM_TYPE pimag = std::sqrt( p.dot(p) );
      Scalar gamma_i = std::sqrt( 1.0 + pimag * pimag );

      // check if curvature photons capable of pair producing can be created.
      if ( !IsProduceCurvPhoton(dt, gamma_i, Rc, rng.uniform()) ) continue;

      Scalar E_ph = SamplePhotonEnergy( gamma_i, Rc );
      Scalar gamma_f = std::max(1.0, gamma_i - E_ph);
      MOM_TYPE pfmag = std::sqrt( gamma_f * gamma_f - 1.0 );

      // create photons
      Vec3<PH_MOM_TYPE> p_ph ( p * (E_ph / pimag) );
      float path_left = _pane.l_coll * f_x( rng.uniform() );
      QNeutral ph_new;
      ph_new.x = x, ph_new.p = p_ph, ph_new.cell = cell;
      ph_new.flag = photon_flag, ph_new.path_left = path_left, ph_new.E = E_ph;
      photons.Append( std::move(ph_new) );
      // only track the event if the parent is tracked
      if ( ParticleTracker::IsTracked(ptc.track_id) ) {
        auto& last_ph = photons.PtcData()[ photons.Number() - 1 ];
        ph_tracker.Track( last_ph );
        auto parentType = particles.Attributes().ptcType;
        pc_tracker.LogThisEvent( timestep, parentType, ptc.track_id, ParticleType::PHOTON, last_ph.track_id );
      }

      num++;

      // primary particle loses energy to gamma rays
      p *= ( pfmag / pimag );
    }
  }

  Logger::print( 0, num, "photons produced in this time step!" );
  InfoCollector::Instance().contents.num_emittedPhotons = num;

}
