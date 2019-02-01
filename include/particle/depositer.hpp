#ifndef _DEPOSITER_HPP_
#define _DEPOSITER_HPP_

namespace knl { enum class shape : unsigned int; }

// TODO check missing 4pi's maybe
namespace particle {
  // NOTE TWJ here may use long double
  template < knl::shape S, typename WJ_field_t, typename ptc_t, typename dq_t, typename grid_t >
  struct depositWJ_t {
    void operator() ( WJ_field_t& WJ, const ptc_t& ptc, const dq_t& dq, const grid_t& grid );
  };

}

#endif
