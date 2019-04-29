#ifndef  _SILO_XX_HPP_
#define  _SILO_XX_HPP_

#include "apt/handle.hpp"
#include "io/silo_operations.hpp"

namespace silo :: traits {
  extern bool display_guard;
  extern int create_mode;
  extern int create_target;
  extern int filetype;
}

typedef struct DBfile DBfile;

namespace silo {
  using DBfileHandle = DBfile*;

  void dbfile_free ( DBfileHandle* p );
  inline DBfileHandle dbfile_null() { return nullptr; };

  enum class Mode : bool { Read = false, Write = true };

  struct file_t : public apt::Handle<DBfileHandle, dbfile_free, dbfile_null>,
                  public SiloPutter<file_t> {};

  template < Mode mode > // use template to prohibit use of Mode::Append in silo::pmpio::open
  file_t open( std::string filename );
}

typedef struct _PMPIO_baton_t PMPIO_baton_t;
namespace mpi { struct Comm; }

namespace silo::pmpio {

  extern int num_files;

  struct DBfileHandle {
    PMPIO_baton_t* baton_h = nullptr;
    DBfile* file_h = nullptr;
  };

  void file_free( DBfileHandle* p );
  inline DBfileHandle file_null() { return {}; }

  struct file_t : public apt::Handle<pmpio::DBfileHandle, pmpio::file_free, pmpio::file_null>,
                  public SiloPutter<pmpio::file_t> {
    // NOTE: This is subtle. This explicit delegated call is needed because otherwise, converting to bool would go as first converting to DBfile*, which is a pointer, and then converting that to bool, which is not what we want. Add explicit in conversion to DBfile* is not preferred because we do want implicit conversion of file_t to DBfile*
    inline operator bool() noexcept {
      return apt::Handle<pmpio::DBfileHandle, pmpio::file_free, pmpio::file_null>::operator bool();
    }

    operator DBfile*() noexcept { return *this ? static_cast<pmpio::DBfileHandle*>(*this)->file_h : nullptr; }

  };

  template < Mode mode >
  file_t open( std::string filename, std::string dirname, const mpi::Comm& comm, int num_files );
}

namespace silo {
  inline void close( file_t& f ) { f.reset(); }
  inline void close( pmpio::file_t& f ) { f.reset(); }
}

#endif
