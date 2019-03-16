#ifndef  _APT_HANDLE_HPP_
#define  _APT_HANDLE_HPP_

#include <memory>

namespace apt {
  template < typename RawHdl,
             void (*Deleter) (RawHdl*),
             RawHdl (*Default) () >
  struct Handle {
  private:
    using shared_ptr = std::shared_ptr<RawHdl>;
    shared_ptr _ptr;
    RawHdl _default = Default();
    RawHdl _fallback = Default();

  public:
    Handle() = default;

    inline void reset() {
      _ptr.reset();
      _fallback = _default;
    }

    inline void reset( RawHdl* raw_handle ) {
      _ptr.reset( raw_handle, Deleter );
    }

    inline void set_fallback_handle( RawHdl fallback ) {
      _fallback = fallback;
    }

    operator RawHdl () const noexcept {
      return use_count() ? *_ptr : _fallback;
    }

    operator RawHdl* () const noexcept {
      return _ptr.get();
    }

    long use_count() const noexcept { return _ptr.use_count(); }
  };
}

#endif
