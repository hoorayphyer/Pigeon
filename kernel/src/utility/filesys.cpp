#include "utility/filesys.hpp"
#include <filesystem>
#include <iostream>

namespace filesystem = std::filesystem;
using std::error_code;

namespace util::fs {
  std::string absolute( std::string dir ) {
    filesystem::path p = dir;
    return filesystem::absolute(p);
  }

  // std::string canonical( std::string dir ) {
  //   filesystem::path p = dir;
  //   return filesystem::canonical(p);
  // }

  std::string& append_slash( std::string& dir ) {
    if ( dir != "" && dir.back() != '/' )
      dir.append("/");
    return dir;
  }

  std::string& remove_slash( std::string& dir ) {
    if ( dir.back() == '/' )
      dir.pop_back();
    return dir;
  }

  void create_directories(std::string dir) {
    // NOTE on some platforms path ending with a slash is not created as a directory
    remove_slash(dir);
    filesystem::path p(dir);
    error_code returnedError;
    if ( !filesystem::create_directories(p, returnedError) && !filesystem::exists(p) ) {
      std::cout << "Directory neither existed nor created: " << dir << std::endl;
    }
    if (returnedError) {
      std::cerr << "ERROR creating directories " << dir << "\nError code = " << returnedError.value() << std::endl;
    }
  }

  bool exists( std::string dir ) {
    filesystem::path p(dir);
    return filesystem::exists( p );
  }

  void remove_all( std::string dir ) {
    filesystem::path p(dir);
    filesystem::remove_all( p );
  }

  void create_directory_symlink( std::string to, std::string new_link ) {
    // NOTE to and new_link must not end with '/'
    remove_slash(to);
    remove_slash(new_link);

    filesystem::path p_to(to);
    filesystem::path p_new_link(new_link);
    filesystem::create_directory_symlink( filesystem::absolute(p_to), filesystem::absolute(p_new_link) );
  }

  void copy_file( std::string target, std::string dest ) {
    filesystem::path target_path(target);
    filesystem::path dest_path(dest);

    error_code returnedError;
    filesystem::copy_file(target_path, dest_path, returnedError);
    if (returnedError) {
      std::cerr << "ERROR copying " << target << " to " << dest << "\nError code = " << returnedError.value() << std::endl;
    }
  }

}
