#include "particle_module/pair_producer.hpp"
#include "utility/rng.hpp"
#include "core/particle.hpp"

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
  template < typename T >
  inline T sample_E_ph(const T& gamma, const T& Rc) {
    return std::min(_pane.E_ph, gamma - 1.0);
  }

  template < typename Iter_el, typename Iter_po,
             typename Tvt, std::size_t DPtc, typename state_t,
             typename Trl = apt::remove_cvref_t<Tvt>,
             typename Ptc = Particle<Tvt, DPtc, state_t> >
  void instant_produce_pairs<Iter_el, Iter_po, Ptc, Trl> ( Iter_el itr_e, Iter_po itr_p, Ptc& ptc, Trl gamma_i, Trl Rc ) {
    {
      // recycle Rc for E_ph
      Rc = sample_E_ph( gamma_i, Rc );
      // primary particle loses energy to gamma rays. ptc.p *= |pf| / |pi|
      ptc.p *= std::sqrt( ( 1.0 - Rc / ( gamma_i - 1.0 ) ) * ( 1.0 - Rc / ( gamma_i + 1.0 ) ) );
    } {
      // recycle Rc for gamma_sec
      Rc /= 2.0;
      // append electron and positron
      auto&& ptc_sec = Particle< Trl, DPtc, apt::remove_cvref_t<state_t> > ( ptc.q, ptc.p * ( std::sqrt( Rc * Rc - 1.0 ) / apt::abs_sq(ptc.p) ) );
      ptc_sec.set<particle::flag::secondary>();
      *(itr_e++) = ptc_sec;
      *(itr_p++) = std::move(ptc_sec);
    }
    // TODO register this pair creation event
    // data.pairCreationEvents.data() [cell] += 1.0;
  }

}

namespace particle {
  template < typename Iter, typename Tvt, std::size_t DPtc, typename state_t,
             typename Trl = apt::remove_cvref_t<Tvt>,
             typename Ptc = Particle<Tvt, DPtc, state_t> >
  void produce_photons<Iter, Ptc, Trl> ( Iter itr_ph, Ptc& ptc, Trl gamma_i, Trl Rc ) {
    // recycle Rc for E_ph
    Rc = sample_E_ph( gamma_i, std::move(Rc) );
    // primary particle loses energy to gamma rays. ptc.p *= |pf| / |pi|
    ptc.p *= std::sqrt( ( 1.0 - Rc / ( gamma_i - 1.0 ) ) * ( 1.0 - Rc / ( gamma_i + 1.0 ) ) );

    *(itr_ph++) = Particle< Trl, DPtc, apt::remove_cvref_t<state_t> >
      ( ptc.q, ptc.p * ( Rc / apt::abs(ptc.p) ) );
  }


  template < typename Iter_el, typename Iter_po,
             typename Tvt, std::size_t DPtc, typename state_t,
             typename Trl = apt::remove_cvref_t<Tvt>,
             typename Ptc = Particle<Tvt, DPtc, state_t> >
  void photon_produce_pairs<Iter_el, Iter_po, Ptc, Trl> ( Iter_el itr_e, Iter_po itr_p, Ptc& photon ) {
    auto&& ptc_sec = Particle< Trl, DPtc, apt::remove_cvref_t<state_t> >
      ( photon.q,
        photon.p * std::sqrt( 0.25 - 1.0 / apt::abs_sq(photon.p) ) );
    ptc_sec.set<particle::flag::secondary>();
    *(itr_e++) = ptc_sec;
    *(itr_p++) = std::move(ptc_sec);

    // TODO register this pair creation event
    // data.pairCreationEvents.data() [cell] += 1.0;

    // erase this photon
    photon.set<particle::flag::empty>();
  }
}
