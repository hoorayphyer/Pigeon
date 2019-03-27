#ifndef  _FIELD_OFFSET_HPP_
#define  _FIELD_OFFSET_HPP_

namespace field {
  struct offset_t {
  private:
    bool _val;
  public:
    constexpr offset_t( bool val = 0 ) noexcept : _val(val) {}
    constexpr offset_t( const offset_t& ) = default;

    constexpr operator double() const noexcept { return static_cast<double>( 0.5 * _val ); }

    constexpr bool operator== ( offset_t other ) const noexcept {
      return _val == other._val;
    }

    constexpr bool operator!= ( offset_t other ) const noexcept {
      return !( *this == other );
    }
  };

}

constexpr field::offset_t INSITU{false}; // right on the grid point
constexpr field::offset_t MIDWAY{true}; // in the middle of an edge connecting two grid points

#endif
