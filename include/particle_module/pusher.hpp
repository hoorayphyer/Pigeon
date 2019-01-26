#ifndef _PARTICLEPUSHER_HPP_
#define _PARTICLEPUSHER_HPP_

#include "vector.hpp"

// TODO check missing 4pi's maybe
namespace particle {
  template < typename T, std::size_t DPtc >
  struct Particle;

  template < typename Tvt, std::size_t DPtc, std::size_t DField,
             typename Trl = apt::remove_cvref_t<Tvt> >
  Vec<Trl,DPtc> update_p( Particle<Tvt,DPtc>& ptc, const Species& sp, Trl dt,
                        const Vec<Trl, DField>& E, const Vec<Trl, DField>& B );

  template < CoordSys CS, typename Tvt, std::size_t DPtc,
             typename Trl = apt::remove_cvref_t<Tvt> >
  Vec<Trl,DPtc> update_q( Particle<Tvt, DPtc>& ptc, const Species& sp, Trl dt );
}




#endif
