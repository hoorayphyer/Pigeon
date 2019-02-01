#ifndef _PAIR_PRODUCER_HPP_
#define _PAIR_PRODUCER_HPP_

class Rng;

namespace particle {
  template < typename ptc_t, typename T >
  bool is_productive_lepton( const ptc_t& ptc, const T& gamma, const T& Rc, Rng& rng ) noexcept;
  template < typename ptc_t, typename T >
  bool is_productive_photon( const ptc_t& photon, const T& B2, Rng& rng ) noexcept;

  template < typename T, typename p_t, typename dp_t >
  T calc_Rc( T dt, const p_t& p, const dp_t& dp ) noexcept;

  template < typename Iter_el, typename Iter_po, typename ptc_t, typename T >
  void instant_produce_pairs( Iter_el itr_e, Iter_po itr_p, ptc_t& ptc, const T& gamma_ptc, T Rc );

  template < typename Iter, typename ptc_t, typename T >
  void produce_photons( Iter itr_photon, ptc_t& ptc, const T& gamma_ptc, T Rc );
  template < typename Iter_el, typename Iter_po, typename ptc_t >
  void photon_produce_pairs( Iter_el itr_e, Iter_po itr_p, ptc_t& photon );
}

#endif
