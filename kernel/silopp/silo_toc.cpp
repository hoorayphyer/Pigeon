#include "silopp/silo_toc.hpp"
#include <silo.h>

namespace silo {
  template < typename file_t >
  std::vector<std::string> Toc<file_t>::toc_dir() {
    std::vector<std::string> res;
    DBtoc* toc = DBGetToc( _dbfile() );
    int num = toc->ndir;
    res.reserve(num);
    for ( int i = 0; i < num; ++i )
      res.emplace_back(toc->dir_names[i]);
    return res;
  }
}

#include "silopp/silo++.hpp"
namespace silo {
  template struct Toc<file_t>;
}
