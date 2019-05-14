#ifndef  _FILESYS_HPP_
#define  _FILESYS_HPP_

#include <string>

namespace fs {
  std::string absolute( std::string dir );

  std::string& append_slash( std::string& dir );

  std::string& remove_slash( std::string& dir );

  void create_directories(std::string dir);

  bool exists( std::string dir );

  void remove_all( std::string dir );

  void create_directory_symlink( std::string to, std::string new_link );

  void copy_file( std::string target, std::string dest );
}

// template wrapper for doing filesystem with mpi. For this to work, include filesys.hpp after mpi. This makes sure no race condition happens during e.g. file creation/deletion
namespace fs {
  template < class Comm, class F >
  void mpido ( const Comm& comm, const F& f ) {
    if ( comm.rank() == 0 ) f();
    comm.barrier();
  }
}

#endif
