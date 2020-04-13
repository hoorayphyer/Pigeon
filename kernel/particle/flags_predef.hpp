#ifndef  _FLAGS_PREDEF_HPP_
#define  _FLAGS_PREDEF_HPP_

namespace particle {
  enum class flag : unsigned int
    { exist = 0, secondary, ignore_force,
      ignore_deposit, ignore_em,
      delimiter, traced };
}

#endif
