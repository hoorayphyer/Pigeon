#ifndef _PAIR_PRODUCER_HPP_
#define _PAIR_PRODUCER_HPP_

class Rng;

namespace particle {
  template < typename Ptc, typename T >
  bool is_productive_lepton( const Ptc& ptc, const T& gamma, const T& Rc, Rng& rng ) noexcept;
  template < typename Ptc, typename T >
  bool is_productive_photon( const Ptc& photon, const T& B2, Rng& rng ) noexcept;

  template < typename Iter_el, typename Iter_po, typename Ptc, typename T >
  void instant_produce_pairs( Iter_el itr_e, Iter_po itr_p, Ptc& ptc, T gamma_ptc, T Rc );

  template < typename Iter, typename Ptc, typename T >
  void produce_photons( Iter itr_photon, Ptc& ptc, T gamma_ptc, T Rc );
  template < typename Iter_el, typename Iter_po, typename Ptc >
  void photon_produce_pairs( Iter_el itr_e, Iter_po itr_p, Ptc& photon );
}

#endif
