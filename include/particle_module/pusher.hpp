#ifndef _PARTICLEPUSHER_HPP_
#define _PARTICLEPUSHER_HPP_

struct Species;
enum class CoordSys;


// TODO check missing 4pi's maybe
namespace particle {
  template < typename Ptc, typename T_field, typename T_dp, typename T >
  T_dp update_p( Ptc& ptc, const Species& sp, const T& dt, const T_field& E, const T_field& B );

  template < CoordSys CS, typename Ptc, typename T_dq, typename T >
  T_dq update_q( Ptc& ptc, const Species& sp, const T& dt );
}




#endif
