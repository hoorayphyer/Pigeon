#ifndef  _SILO_XX_HPP_
#define  _SILO_XX_HPP_

#include "apt/handle.hpp"
#include "silopp/silo_operations.hpp"

typedef struct DBfile DBfile;

namespace silo {
  using DBfileHandle = DBfile*;

  void dbfile_free ( DBfileHandle* p );
  inline DBfileHandle dbfile_null() { return nullptr; };

  enum class Mode : bool { Read = false, Write = true };

  struct file_t : public apt::Handle<DBfileHandle, dbfile_free, dbfile_null>,
                  public SiloPutter<file_t> {};

  file_t open( std::string filename, Mode mode );
  inline void close( file_t& f ) { f.reset(); }
}

#endif
