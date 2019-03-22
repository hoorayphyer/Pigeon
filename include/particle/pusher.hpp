#ifndef _PARTICLE_PUSHER_HPP_
#define _PARTICLE_PUSHER_HPP_

#include "kernel/coordsys_predef.hpp"
#include "particle/particle_expression.hpp"
#include "apt/vec_expression.hpp"

namespace particle {
  template < typename dp_t, typename Ptc, typename Field, typename T >
  dp_t update_p( PtcExpression<Ptc>& ptc, T dt, unsigned int mass_x,
                 const apt::VecExpression<Field>& E, const apt::VecExpression<Field>& B );

  template < species sp, knl::coordsys CS, typename dq_t, typename Ptc, typename T >
  dq_t update_q( PtcExpression<Ptc>& ptc, T dt, bool is_massive );
}




#endif
