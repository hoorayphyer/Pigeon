#ifndef _PAIR_PRODUCER_HPP_
#define _PAIR_PRODUCER_HPP_

#include "particle/particle_expression.hpp"
#include "apt/vec_expression.hpp"

namespace util { template < typename > struct Rng; }

namespace particle {
  template < typename Ptc, typename T = typename Ptc::value_type >
  bool is_productive_lepton( const PtcExpression<Ptc>& ptc, const T& gamma,
                             const T& Rc, util::Rng<T>& rng ) noexcept;
  template < typename Ptc, typename T = typename Ptc::value_type >
  bool is_productive_photon( const PtcExpression<Ptc>& photon, const T& B2,
                             util::Rng<T>& rng ) noexcept;

  template < typename T, typename p_t, typename dp_t >
  T calc_Rc( T dt, const apt::VecExpression<p_t>& p, const apt::VecExpression<dp_t>& dp ) noexcept;

  template < typename BackInsertIter, typename Ptc, typename T = typename Ptc::value_type >
  void instant_produce_pairs( BackInsertIter itr_e, BackInsertIter itr_p,
                              PtcExpression<Ptc>& ptc, const T& gamma_ptc, T Rc );
  template < typename BackInsertIter, typename Ptc, typename T = typename Ptc::value_type >
  void produce_photons( BackInsertIter itr_photon, PtcExpression<Ptc>& ptc,
                        const T& gamma_ptc, T Rc );
  template < typename BackInsertIter, typename Ptc >
  void photon_produce_pairs( BackInsertIter itr_e, BackInsertIter itr_p, PtcExpression<Ptc>& photon );
}

#endif
