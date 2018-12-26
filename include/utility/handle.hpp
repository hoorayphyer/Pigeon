#ifndef  _HANDLE_HPP_
#define  _HANDLE_HPP_
#include <memory>

// TODOL this struct is not very type safe, i.e. cannot trigger type error at compile time
struct Handle {
private:
  std::shared_ptr<void> _ptr = nullptr; // type erasure
public:
  Handle() = default;

  template < typename RawHdl, typename Deleter >
  Handle( RawHdl* raw_handle, Deleter d ) : _ptr( raw_handle, d ) {}

  inline explicit operator bool () const noexcept {
    return static_cast<bool>(_ptr);
  }

  template < typename RawHdl >
  operator  RawHdl () const {
    return *std::static_pointer_cast<RawHdl>(_ptr);
  }

  template < typename RawHdl >
  operator RawHdl& () {
    return *std::static_pointer_cast<RawHdl>(_ptr);
  }

  inline void reset() { _ptr.reset(); }

  template < typename RawHdl >
  inline void reset( RawHdl* raw_handle ) { _ptr.reset(raw_handle); } // TODOL deleter should remain the same

  bool operator==( Handle other ) const {
    return _ptr.get() == other._ptr.get(); // TODOL check correctness and how to check type equality?
  }

  inline bool operator!=( Handle other ) const { return !operator==(other); }

};

const Handle nullhdl; // to mimick nullptr // TODOL use some nullptr_t kind of std template?

#endif
