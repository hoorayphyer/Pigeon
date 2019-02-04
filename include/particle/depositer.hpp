#ifndef _DEPOSITER_HPP_
#define _DEPOSITER_HPP_

namespace knl { enum class shape : unsigned int; }

// TODO check missing 4pi's maybe
namespace particle {
  // NOTE TWJ here may use long double
  template < knl::shape S, typename Field, typename Ptc, typename Vec, typename Grid >
  void depositWJ ( Field& WJ, const Ptc& ptc, const Vec& dq, const Grid& grid );
}

#endif
