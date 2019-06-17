#ifndef _TIMER_HPP_
#define _TIMER_HPP_
#include <any>
#include <string>

namespace tmr {
  struct Timestamp {
  private:
    std::any _t;
    std::string _unit{"ms"}; // possible values are s, ms, us, ns
  public:
    Timestamp();

    double lapse(std::string unit) const;
    inline double lapse() const { return lapse(_unit); }

    inline void set_unit(std::string unit) noexcept {
      _unit = std::move(unit);
      if ( _unit != "ms" && _unit != "us" && _unit != "ns" && _unit != "s" ) _unit = "ms";
    }
    inline const std::string& unit() const noexcept { return _unit; }

    double operator- ( const Timestamp& other );
  };

}

#endif
