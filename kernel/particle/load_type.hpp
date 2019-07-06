#ifndef _PARTICLE_LOAD_TYPE_HPP_
#define _PARTICLE_LOAD_TYPE_HPP_

namespace particle {
  // use long long instead of unsigned long long to avoid sutble bug like
  // int a = 1; unsigned long b = 100; auto x = a - b;
  // in which up conversion will do auto x = (unsigned long) a - b, resulting in x being (unsigned long)(-99)
  using load_t = long long;
}

#endif
