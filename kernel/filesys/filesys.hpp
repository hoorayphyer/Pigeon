#ifndef  _FILESYS_HPP_
#define  _FILESYS_HPP_

#include <string>
#include <any>

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

namespace fs {
  struct directory_iterator {
  private:
    std::any _itr;

  public:
    using value_type = std::string;
    using difference_type =	std::ptrdiff_t;
    using pointer =	const std::string*;
    using reference =	const std::string&;
    using iterator_category =	std::input_iterator_tag;

    directory_iterator();
    directory_iterator( std::string );
    directory_iterator(const directory_iterator&) = default;
    directory_iterator(directory_iterator&&) = default;

    directory_iterator& operator=( const directory_iterator& ) = default;
    directory_iterator& operator=( directory_iterator&& ) = default;

    bool operator!=( const directory_iterator& ) const;

    value_type operator*() const;
    directory_iterator& operator++();

    directory_iterator begin() const noexcept { return *this; }
    directory_iterator end() const noexcept;
  };
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
