#ifndef _APT_FOREACH_HPP_
#define _APT_FOREACH_HPP_

#include <utility> // for std::forward

namespace apt {
  template < int Begin, int End, typename Func, typename... Args >
  constexpr void foreach( const Func& f, Args&&... args  ) noexcept {
    static_assert( Begin <= End );
    // TODO add bounds check
    // static_assert( (... && ( End - Begin <= Args::NDim ) ) );
    for ( auto i = Begin; i < End; ++i )
      f( std::forward<Args>(args)[i]... );
  }

}

#endif
