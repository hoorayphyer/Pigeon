#include "particle/pair_producer.hpp"
#include "utility/rng.hpp"
#include "particle/particle.hpp"
#include "apt/numeric.hpp"

namespace particle {
  template < typename Ptc, typename T >
  bool is_productive_lepton( const Ptc& ptc, const T& gamma, const T& Rc, Rng& rng ) noexcept {
    return
      _pane.in_productive_zone( ptc.q ) &&
      gamma > _pane.gamma_off &&
      gamma > _pane.K_thr *  std::cbrt(Rc) &&
      // prob_fiducial = K_curv_em_rate * dt
      rng.uniform() < _pane.prob_fiducial * gamma  / Rc;
  }
}

namespace particle {
  namespace opacity {
    // prob_mag_conv = dt / mfp_mag_conv
    template < typename T >
    inline bool mag_conv( const T& B2, T randnum ) noexcept {
      return B2 > _pane.B_magconv * _pane.B_magconv ? ( randnum < prob_mag_conv ) : false;
    }

    // template < typename T >
    // inline T f_x ( T x ) {
    //   // distribution of x*exp(-x^2/2), which peaks at x = 1.
    //   return std::sqrt( -2.0 * std::log(x) );
    // }

    // prob_mag_conv = dt / mfp_ph_ph
    template < typename T >
    inline bool ph_ph( T randnum ) noexcept {
      // TODO double check this implementation, it is not equivalent because the original one has some sort of gaussian in it. Use Monte Carlo
      return randnum < prob_ph_ph;
    }
  }

  template < typename Ptc, typename T >
  bool is_productive_photon( const Ptc& photon, const T& B2, Rng& rng ) noexcept {
    return opacity::magnetic_convert(B2, rng.uniform() ) || opacity::ph_ph( rng.uniform() );
  }
}

namespace particle {
  template < typename Tvt, std::size_t DPtc,
             typename Trl = apt::remove_cvref_t<Tvt>,
             typename T_p = Vec<Tvt, DPtc>,
             typename T_dp = Vec<Trl, DPtc>
             >
  Trl calc_Rc<Trl, T_p, T_dp>( Trl dt, const T_p& p, const T_dp& dp ) noexcept {
    // TODO don't use a uniform number for Rc
    return 1.0;
  }

  // float ParticlePusher::CalculateRc( Scalar dt, const Vec3<MOM_TYPE> &p, const Vec3<MOM_TYPE> &dp) const {
  //   // find momentum at half time step
  //   Vec3<MOM_TYPE> phalf( p + dp * 0.5 );
  //   Vec3<MOM_TYPE> v( phalf / std::sqrt( 1.0 + phalf.dot(phalf) ) );
  //   Scalar vv = v.dot( v );
  //   Vec3<MOM_TYPE> a( dp / dt ); // a is for now force, will be converted to dv/dt
  //   // convert a to dv/dt
  //   a = ( a - v * ( v.dot(a) ) ) * std::sqrt( 1.0 - vv );
  //   Scalar va = v.dot( a ); // get the real v dot a
  //   return vv / std::max( std::sqrt( a.dot(a) - va * va / vv ), 1e-6 ); // in case denominator becomes zero
  // }

  // float ParticlePusher::GetDipolarRc(const Scalar &r_sph, const Scalar &cos_th, const Scalar &phi) const {
  //   Scalar sin_th = std::sqrt( 1.0 - cos_th * cos_th );
  //   Scalar tmp1 = 1.0 + cos_th * cos_th;
  //   Scalar tmp2 = 3.0 * tmp1 - 2.0;
  //   return r_sph * tmp2 * std::sqrt(tmp2) / ( 3.0 * tmp1 * sin_th );
  // }

}

namespace particle {
  template < typename T >
  inline T sample_E_ph(const T& gamma, const T& Rc) {
    return std::min(_pane.E_ph, gamma - 1.0);
  }

  template < typename Iter_el, typename Iter_po,
             typename Tvt, std::size_t DPtc,
             typename Trl = apt::remove_cvref_t<Tvt>,
             typename Ptc = Particle<Tvt, DPtc> >
  void instant_produce_pairs<Iter_el, Iter_po, Ptc, Trl> ( Iter_el itr_e, Iter_po itr_p, Ptc& ptc, const Trl& gamma_i, Trl Rc ) {
    {
      // recycle Rc for E_ph
      Rc = sample_E_ph( gamma_i, Rc );
      // primary particle loses energy to gamma rays. ptc.p *= |pf| / |pi|
      ptc.p *= std::sqrt( ( 1.0 - Rc / ( gamma_i - 1.0 ) ) * ( 1.0 - Rc / ( gamma_i + 1.0 ) ) );
    } {
      // recycle Rc for gamma_sec
      Rc /= 2.0;
      // append electron and positron
      auto&& ptc_sec = Particle< Trl, DPtc > ( ptc.q,
                                               ptc.p * ( std::sqrt( Rc * Rc - 1.0 ) / apt::abs_sq(ptc.p) ),
                                               flag::secondary);
      ptc_sec.set(species::electron);
      *(itr_e++) = ptc_sec;
      ptc_sec.set(species::positron);
      *(itr_p++) = std::move(ptc_sec);
    }
    // TODO register this pair creation event
    // data.pairCreationEvents.data() [cell] += 1.0;
  }

}

namespace particle {
  template < typename Iter, typename Tvt, std::size_t DPtc,
             typename Trl = apt::remove_cvref_t<Tvt>,
             typename Ptc = Particle<Tvt, DPtc> >
  void produce_photons<Iter, Ptc, Trl> ( Iter itr_ph, Ptc& ptc, const Trl& gamma_i, Trl Rc ) {
    // recycle Rc for E_ph
    Rc = sample_E_ph( gamma_i, std::move(Rc) );
    // primary particle loses energy to gamma rays. ptc.p *= |pf| / |pi|
    ptc.p *= std::sqrt( ( 1.0 - Rc / ( gamma_i - 1.0 ) ) * ( 1.0 - Rc / ( gamma_i + 1.0 ) ) );

    *(itr_ph++) = Particle< Trl, DPtc > ( ptc.q,
                                          ptc.p * ( Rc / apt::abs(ptc.p) ),
                                          species::photon );
  }


  template < typename Iter_el, typename Iter_po,
             typename Tvt, std::size_t DPtc,
             typename Trl = apt::remove_cvref_t<Tvt>,
             typename Ptc = Particle<Tvt, DPtc> >
  void photon_produce_pairs<Iter_el, Iter_po, Ptc, Trl> ( Iter_el itr_e, Iter_po itr_p, Ptc& photon ) {
    auto&& ptc_sec = Particle< Trl, DPtc >
      ( photon.q,
        photon.p * std::sqrt( 0.25 - 1.0 / apt::abs_sq(photon.p) ),
        flag::secondary );
    ptc_sec.set(species::electron);
    *(itr_e++) = ptc_sec;
    ptc_sec.set(species::positron);
    *(itr_p++) = std::move(ptc_sec);

    // TODO register this pair creation event
    // data.pairCreationEvents.data() [cell] += 1.0;

    // void this photon
    photon.set(flag::empty);
  }
}
