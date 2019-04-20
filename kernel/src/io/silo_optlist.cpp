#include "io/silo_optlist.hpp"
#include <silo.h>

namespace silo {
  constexpr bool float_opt( int option ) noexcept {
    switch( option ) {
    case DBOPT_TIME: return true;
    default: return false;
    }
  }

  constexpr bool double_opt( int option ) noexcept {
    switch( option ) {
    case DBOPT_DTIME: return true;
    default: return false;
    }
  }

  constexpr bool int_opt( int option ) noexcept {
    switch( option ) {
    case DBOPT_CYCLE: return true;
    default: return false;
    }
  }

  OptVal::OptVal( int id ) {
    if ( int_opt(id) )  val = (int)0;
    else if ( float_opt(id) ) val = (float)0;
    else if ( double_opt(id) ) val = (double)0;
  }
};

namespace silo {
  void optlist_free (DBoptlist * optlist) {
    DBFreeOptlist(optlist);
  }

  OptList::operator DBoptlist* () {
    _p.reset(DBMakeOptlist( size() ));

    for ( const auto& option : *this ) {
      auto id = option.first;
      const auto& val = option.second.val;
      std::visit( [id, ptr = _p.get()](auto& v)
                  { DBAddOption(ptr, id, (void*)&v); }, val );
    }

    return _p.get();
  }
}
