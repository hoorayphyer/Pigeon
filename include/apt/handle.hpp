#ifndef  _APT_HANDLE_HPP_
#define  _APT_HANDLE_HPP_

#include <memory>

namespace apt {
  template < typename RawHdl, void (*Deleter) ( RawHdl* ) >
  struct Handle : public std::shared_ptr<RawHdl> {
    Handle() = default;

    Handle( RawHdl* raw_handle ) {
      reset( raw_handle );
    }

    inline void reset() {
      std::shared_ptr<RawHdl>::reset();
    }

    inline void reset( RawHdl* raw_handle ) {
      std::shared_ptr<RawHdl>::reset( raw_handle, Deleter );
    }

    operator RawHdl () const noexcept { return *this; }

    operator RawHdl* () const noexcept { return this->get(); }
  };
}

#endif
