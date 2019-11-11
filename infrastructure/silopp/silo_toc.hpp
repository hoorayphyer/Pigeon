#ifndef _SILO_TOC_
#define _SILO_TOC_

#include <vector>
#include <string>

typedef struct DBfile DBfile;

// NOTE DBSetDir will invalidate DBtoc* obtained by calling DBGetToc()
namespace silo {
  template < typename file_t >
  struct Toc {
  private:
    inline DBfile* _dbfile() noexcept { return static_cast<file_t&>(*this).operator DBfile* (); }

    std::vector<std::string> toc_impl( int num, char** names );

  public:
    std::vector<std::string> toc_dir();
    // std::vector<std::string> toc_array();
    std::vector<std::string> toc_var();

  };

}

#endif
