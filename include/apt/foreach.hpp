#ifndef _APT_FOREACH_HPP_
#define _APT_FOREACH_HPP_

#include <utility> // for std::forward

namespace apt {
  template < std::size_t Begin, std::size_t End, typename Func, typename... Args >
  constexpr void foreach( const Func& f, Args&&... args  ) noexcept {
    static_assert( Begin <= End );
    if constexpr ( Begin == End ) return;
    else {
      f( std::get<Begin>(std::forward<Args>(args))... );
      return foreach<Begin+1, End>( f, std::forward<Args>(args)... );
    }
  }

}

#endif
