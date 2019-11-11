#ifndef  _SILO_NAVIGATION_HPP_
#define  _SILO_NAVIGATION_HPP_
#include <string>

typedef struct DBfile DBfile;

namespace silo {

  template < typename file_t >
  struct Navigator {
  private:
    inline DBfile* _dbfile() noexcept { return static_cast<file_t&>(*this).operator DBfile* (); }
  public:
    void mkdir( const std::string& path );
    void cd ( const std::string& path );
    bool exists( const std::string& path );

    // cd to the path and if path is not existent, create it first
    void mkcd( const std::string& path );
  };

}

#endif
