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
    void mkdir( std::string path );
    void cd ( std::string path );
    bool exists( std::string path );
  };

}

#endif
