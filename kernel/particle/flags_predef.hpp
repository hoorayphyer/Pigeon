#ifndef  _FLAGS_PREDEF_HPP_
#define  _FLAGS_PREDEF_HPP_

namespace particle {
  // 16 flags
  enum class flag : unsigned int
    { exist = 0,
      secondary,
      traced,
      ignore_force,
      ignore_deposit,
      _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15
    };

  constexpr bool is_reserved_flag( flag f ) noexcept {
    return static_cast<unsigned int>(f) < 5;
  }
}

#endif
