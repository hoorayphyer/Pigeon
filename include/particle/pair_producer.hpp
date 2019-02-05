#ifndef _PAIR_PRODUCER_HPP_
#define _PAIR_PRODUCER_HPP_

namespace util {
  template < typename > struct Rng; }


namespace particle {
  template < typename Ptc, typename T >
  bool is_productive_lepton( const Ptc& ptc, const T& gamma,
                             const T& Rc, util::Rng<T>& rng ) noexcept;
  template < typename Ptc, typename T >
  bool is_productive_photon( const Ptc& photon, const T& B2,
                             util::Rng<T>& rng ) noexcept;

  template < typename T, typename Vec_p_t, typename Vec_dp_t >
  T calc_Rc( T dt, const Vec_p_t& p, const Vec_dp_t& dp ) noexcept;

  template < typename Iter, typename Ptc, typename T >
  void instant_produce_pairs( Iter itr_e, Iter itr_p, Ptc& ptc, const T& gamma_ptc, T Rc );

  template < typename Iter, typename Ptc, typename T >
  void produce_photons( Iter itr_photon, Ptc& ptc, const T& gamma_ptc, T Rc );
  template < typename Iter, typename Ptc >
  void photon_produce_pairs( Iter itr_e, Iter itr_p, Ptc& photon );
}

#endif
