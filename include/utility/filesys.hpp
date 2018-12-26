#ifndef  _FILESYS_HPP_
#define  _FILESYS_HPP_

#include <string>

namespace filesys {
  std::string append_slash( std::string dir );

  void create_directories(std::string dir);

  bool exists( std::string dir );

  void remove_all( std::string dir );

  void create_directory_symlink( std::string to, std::string new_link );

  void copy_file( std::string target, std::string dest );
}

#endif
