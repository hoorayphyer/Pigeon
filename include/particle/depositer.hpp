#ifndef _DEPOSITER_HPP_
#define _DEPOSITER_HPP_

// TODO check missing 4pi's maybe
namespace particle {
  // NOTE TWJ here may use long double
  template < sf::shape S, typename T_WJ_Field, typename Ptc, typename T_dq, typename T_Grid >
  void depositWJ( T_WJ_Field& WJ, const Ptc& ptc, const T_dq& dq, const T_Grid& grid );
}

#endif
