#ifndef _PARTICLEPUSHER_HPP_
#define _PARTICLEPUSHER_HPP_

namespace particle { enum class species : unsigned char; }
namespace knl { enum class coordsys_t : unsigned char; }

// TODO check missing 4pi's maybe
namespace particle {
  template < species sp, typename Ptc, typename Field, typename Vec, typename T >
  Vec update_p( Ptc& ptc, const T& dt, const Field& E, const Field& B );

  template < species sp, knl::coordsys_t CS, typename Ptc, typename Vec, typename T >
  Vec update_q( Ptc& ptc, const T& dt );
}




#endif
