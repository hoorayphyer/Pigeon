#include "utility/filesys.hpp"
#include <filesystem>
#include <iostream>

namespace filesys {
  std::string append_slash( std::string dir ) {
    if ( dir != "" && dir.back() != '/' )
      dir.append('/');
    return dir;
  }

  void create_directories(std::string dir) {
    boost::filesystem::path p(dir);
    boost::system::error_code returnedError;
    boost::filesystem::create_directories(p, returnedError);
    if (returnedError) {
      std::cerr << "ERROR creating directories " << dir << "\nError code = " << returnedError.value() << std::endl;
    }
  }

  bool exists( std::string dir ) {
    boost::filesystem::path p(dir);
    return boost::filesystem::exists( p );
  }

  void remove_all( std::string dir ) {
    boost::filesystem::path p(dir);
    boost::filesystem::remove_all( p );
  }

  void create_directory_symlink( std::string to, std::string new_link ) {
    // NOTE to and new_link must not end with '/'
    if (to.back() == '/')
      to.pop_back();

    if (new_link.back() == '/')
      new_link.pop_back();

    boost::filesystem::path p_to(to);
    boost::filesystem::path p_new_link(new_link);
    boost::filesystem::create_directory_symlink( boost::filesystem::absolute(p_to), boost::filesystem::absolute(p_new_link) );

  }

  void copy_file( std::string target, std::string dest ) {
    boost::filesystem::path target_path(target);
    boost::filesystem::path dest_path(dest);

    boost::system::error_code returnedError;
    boost::filesystem::copy_file(target_path, dest_path, returnedError);
    if (returnedError) {
      std::cerr << "ERROR copying " << target << " to " << dest << "\nError code = " << returnedError.value() << std::endl;
    }
  }

}
