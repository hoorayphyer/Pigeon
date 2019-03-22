#include "particle/pair_producer.cpp"
#include "traits.hpp"
#include "particle/virtual_particle.hpp"
#include "particle/array.hpp"

using namespace traits;

namespace particle {
  using vPtc = vParticle< real_t, DPtc, ptc_state_t >;
  using BackInsertIter = std::back_insert_iterator< array<real_t,DPtc,ptc_state_t>>;

  template bool
  is_productive_lepton<vPtc, real_t>( const PtcExpression<vPtc>& ptc,
                                      const real_t& gamma,
                                      const real_t& Rc,
                                      util::Rng<real_t>& rng ) noexcept;

  template bool
  is_productive_photon<vPtc, real_t>( const PtcExpression<vPtc>& photon,
                                      const real_t& B2,
                                      util::Rng<real_t>& rng ) noexcept;

  template real_t
  calc_Rc<real_t,
          apt::vVec<real_t,DPtc>,
          apt::Vec<real_t,DPtc>> ( real_t dt,
                                   const apt::VecExpression<apt::vVec<real_t,DPtc>>& p,
                                   const apt::VecExpression<apt::Vec<real_t,DPtc>>& dp ) noexcept;

  template void
  instant_produce_pairs< BackInsertIter,
                         vPtc,
                         real_t> ( BackInsertIter itr_e,
                                   BackInsertIter itr_p,
                                   PtcExpression<vPtc>& ptc,
                                   const real_t& gamma_ptc,
                                   real_t Rc );

  template void
  produce_photons< BackInsertIter,
                   vPtc,
                   real_t > ( BackInsertIter itr_photon,
                              PtcExpression<vPtc>& ptc,
                              const real_t& gamma_ptc,
                              real_t Rc );

  template void
  photon_produce_pairs<BackInsertIter,
                       vPtc> ( BackInsertIter itr_e,
                               BackInsertIter itr_p,
                               PtcExpression<vPtc>& photon );
}
