#ifndef _TIMER_HPP_
#define _TIMER_HPP_
#include <any>

namespace tmr {
  struct TimeStamp {
  private:
    std::any _t;
  public:
    TimeStamp();
    double lapse();
    double operator- ( const TimeStamp& other );
  };
}

#endif
