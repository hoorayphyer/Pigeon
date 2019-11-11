#ifndef  _SILO_XX_HPP_
#define  _SILO_XX_HPP_

#include "apt/handle.hpp"
#include "silopp/silo_operations.hpp"
#include "silopp/silo_navigation.hpp"
#include "silopp/silo_reader.hpp"
#include "silopp/silo_toc.hpp"

typedef struct DBfile DBfile;

namespace silo {
  using DBfileHandle = DBfile*;

  void dbfile_free ( DBfileHandle* p );
  inline DBfileHandle dbfile_null() { return nullptr; };

  enum class Mode : bool { Read = false, Write = true };

  struct file_t : public apt::Handle<DBfileHandle, dbfile_free, dbfile_null>,
                  public Operations<file_t>,
                  public Reader<file_t>,
                  public Navigator<file_t>,
                  public Toc<file_t> {};

  file_t open( std::string filename, Mode mode );
  inline void close( file_t& f ) { f.reset(); }
}

namespace silo {
  std::string errmsg();
}

#endif
