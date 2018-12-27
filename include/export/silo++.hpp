#ifndef  _SILO_XX_HPP_
#define  _SILO_XX_HPP_

#include "utility/handle.hpp"
#include <string>
#include <memory>

namespace mpi {
  struct Comm;
}

namespace silo :: traits {
  extern bool display_guard;
  extern int create_mode;
  extern int create_target;
  extern int filetype;
}

namespace silo {
  enum class Mode : char { Read = 0, Write, Append };

  struct DBfile_t {
    Handle _baton;
    Handle file_h;
    ~DBfile_t();
  };

  void close( DBfile_t );

  template < Mode mode >
  DBfile_t open( std::string filename );

  // extern std::function<void(DBfile*)> put_mesh;
  // extern std::function<void(DBfile*)> put_var;
  // extern std::function<void(DBfile*, int timestep)> put_multimesh;
  // extern std::function<void(DBfile*, int timestep, std::string name_sth, const int type)> put_multivar;


}

namespace silo::pmpio {
  extern int num_files;

  template < Mode mode >
  DBfile_t open( std::string dirname, const mpi::Comm& comm );
}

#endif
