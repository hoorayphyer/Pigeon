#ifndef  _SILO_OPERATIONS_HPP_
#define  _SILO_OPERATIONS_HPP_
#include <string>
#include <vector>
#include "io/silo_optlist.hpp"

typedef struct DBfile DBfile;

namespace silo {

  template < typename file_t >
  struct SiloPutter {
  private:
    DBfile* _dbfile() noexcept { return static_cast<file_t&>(*this)._dbfile(); }

  public:
    template < typename T >
    void put_mesh( std::string meshname, const std::vector<std::vector<T>>& coords, const OptList& optlist );

    template < typename T >
    void put_var( std::string varname, std::string meshname, const T* vardata, const std::vector<int>& dims );

    void put_multimesh( std::string multimeshname, int nblock, std::string file_ns, std::string block_ns, OptList optlist );

    void put_multivar( std::string multivarname, int nblock, std::string file_ns, std::string block_ns, OptList optlist );
  };
}

#endif
