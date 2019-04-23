#include "io/silo_optlist.hpp"
#include <silo.h>

namespace silo {
  constexpr bool float_opt( int option ) noexcept {
    switch( option ) {
    case DBOPT_TIME:
      return true;
    default:
      return false;
    }
  }

  constexpr bool double_opt( int option ) noexcept {
    switch( option ) {
    case DBOPT_DTIME:
      return true;
    default:
      return false;
    }
  }

  constexpr bool int_opt( int option ) noexcept {
    switch( option ) {
    case DBOPT_CYCLE:
      return true;
    default:
      return false;
    }
  }

  constexpr bool int_array_opt( int option ) noexcept {
    switch( option ) {
    case DBOPT_LO_OFFSET:
    case DBOPT_HI_OFFSET:
      return true;
    default:
      return false;
    }
  }

  OptVal::OptVal( int id ) {
    if ( int_opt(id) )  _numeric = (int)0;
    else if ( float_opt(id) ) _numeric = (float)0;
    else if ( double_opt(id) ) _numeric = (double)0;
    else if ( int_array_opt(id) ) _array = std::vector<int>{};
  }
};

namespace silo {
  void optlist_free (DBoptlist * optlist) {
    DBFreeOptlist(optlist);
  }

  OptList::operator DBoptlist* () {
    _p.reset(DBMakeOptlist( size() ));

    for ( auto& [id, val] : *this ) {
      if ( int_array_opt(id) ) {
        DBAddOption(_p.get(), id, val._array.data());
      } else {
        std::visit( [id, ptr = _p.get()](auto& v)
                    { DBAddOption(ptr, id, (void*)&v); }, val._numeric );
      }
    }

    return _p.get();
  }
}
