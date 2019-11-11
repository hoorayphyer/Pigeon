#include "silopp/silo_navigation.hpp"
#include <silo.h>

namespace silo {
  template < typename file_t >
  void Navigator<file_t>::mkdir( const std::string& path ) {
    DBMkDir(_dbfile(), path.c_str());
  }

  template < typename file_t >
  void Navigator<file_t>::cd( const std::string& path ) {
    DBSetDir(_dbfile(), path.c_str());
  }

  template < typename file_t >
  bool Navigator<file_t>::exists( const std::string& path ) {
    return DBInqVarExists(_dbfile(), path.c_str());
  }

  template < typename file_t >
  void Navigator<file_t>::mkcd( const std::string& path ) {
    if ( !exists(path) ) mkdir(path);
    cd(path);
  }
}


#include "silopp/silo++.hpp"
namespace silo {
  template struct Navigator<file_t>;
}
