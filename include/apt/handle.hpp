#ifndef  _APT_HANDLE_HPP_
#define  _APT_HANDLE_HPP_

#include <memory>

namespace apt {
  template < typename RawHdl, void (*Deleter) ( RawHdl* ) >
  struct Handle : public std::shared_ptr<RawHdl> {
    Handle() : std::shared_ptr<RawHdl>( nullptr ) {}

    Handle( RawHdl* raw_handle )
      : std::shared_ptr<RawHdl>( raw_handle, Deleter ) {}

    operator RawHdl () const noexcept { return *this; }

    operator RawHdl* () const noexcept { return this->get(); }
  };
}

#endif
