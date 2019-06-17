#ifndef _TIMER_HPP_
#define _TIMER_HPP_
#include <any>
#include <string>

namespace tmr {
  using TDur = double;
  struct Duration {
  private:
    std::string _unit{"ms"}; // possible values are s, ms, us, ns
    const TDur _val{};

    Duration( TDur val ) noexcept : _val(std::move(val)) {}

  public:
    friend struct Timestamp;

    inline Duration& in_units_of(std::string unit) noexcept {
      _unit = std::move(unit);
      if ( _unit != "ms" && _unit != "us" && _unit != "ns" && _unit != "s" ) _unit = "ms";
      return *this;
    }

    inline const std::string& unit() const noexcept { return _unit; }
    TDur val(std::string unit) const noexcept;
    inline TDur val() const noexcept { return val(_unit); }
  };

  template < typename OStream >
  OStream& operator<< ( OStream& os, const Duration& dur ) {
    os << dur.val(dur.unit()) << " " << dur.unit();
    return os;
  }

  struct Timestamp {
  private:
    std::any _t;
  public:
    Timestamp();
    Duration lapse() const;
    Duration operator- ( const Timestamp& other );
  };

}

#endif
