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

  public:
    Handle() = default;

    inline void reset() { _ptr.reset(); }

    inline void reset( RawHdl* raw_handle ) {
      _ptr.reset( raw_handle, Deleter );
    }

    operator RawHdl () const noexcept {
      return use_count() ? *_ptr : Default();
    }

    operator RawHdl* () const noexcept {
      return _ptr.get();
    }

    long use_count() const noexcept { return _ptr.use_count(); }
  };
}

#endif
