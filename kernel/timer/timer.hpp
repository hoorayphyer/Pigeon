#ifndef _TIMER_HPP_
#define _TIMER_HPP_
#include <any>

namespace tmr {
  struct Timestamp {
  private:
    std::any _t;
  public:
    Timestamp();
    double lapse();
    double operator- ( const Timestamp& other );
  };

}

#endif
