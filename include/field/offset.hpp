#ifndef  _FIELD_OFFSET_HPP_
#define  _FIELD_OFFSET_HPP_

namespace field {

  struct offset_t {
  private:
    bool _val;
  public:
    constexpr offset_t( bool val ) noexcept : _val(val) {}
    constexpr offset_t( const offset_t& ) = default;

    constexpr operator bool() noexcept { return _val; }

    template < typename T >
    constexpr operator T() noexcept { return static_cast<T>( 0.5 * _val ); }
  };

  constexpr offset_t INSITU{false}; // right on the grid point
  constexpr offset_t MIDWAY{true}; // in the middle of an edge connecting two grid points
}

#endif
