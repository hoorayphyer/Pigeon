#ifndef  _SILO_OPERATIONS_HPP_
#define  _SILO_OPERATIONS_HPP_
#include <string>
#include <vector>
#include "silopp/silo_optlist.hpp"

typedef struct DBfile DBfile;

namespace silo {
  enum class MeshType : char { Rect=0, Curv };
}

namespace silo {

  template < typename file_t >
  struct SiloPutter {
  private:
    inline DBfile* _dbfile() noexcept { return static_cast<file_t&>(*this).operator DBfile* (); }

  public:
    template < typename T >
    void put_mesh( std::string meshname, const std::vector<std::vector<T>>& coords, MeshType mt, OptList optlist = {} );

    // put scalar field
    template < typename T >
    void put_var( std::string varname, std::string meshname, const T* vardata, const std::vector<int>& dims, OptList optlist = {} );

    // put vector or tensor field
    template < typename T >
    void put_var( std::string varname, std::string meshname, const std::vector<const T*>& vardata, const std::vector<int>& dims, OptList optlist = {} );

    void put_multimesh( std::string multimeshname, int nblock, std::string file_ns, std::string block_ns, MeshType mt, OptList optlist = {} );

    void put_multivar( std::string multivarname, int nblock, std::string file_ns, std::string block_ns, OptList optlist = {} );
  };
}

#endif
