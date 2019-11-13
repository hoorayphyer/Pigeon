#include "silopp/silo++.hpp"
#include "silo.h"

namespace silo :: traits {
  bool display_guard = true;
  int create_mode = DB_CLOBBER;
  int create_target = DB_LOCAL;
  int filetype = DB_HDF5;

  // TODO these are global settings
  // DBSetFriendlyHDF5Names(1);
  // DBSetCompression("METHOD=GZIP");
}

namespace silo {
  void dbfile_free ( DBfileHandle* p ) {
    if (p && *p) DBClose(*p);
    delete p;
  }

  file_t open(  std::string filename, Mode mode ) {
    file_t dbfile;

    if  ( Mode::Read == mode )
      dbfile.reset( new DBfileHandle(DBOpen( filename.c_str(), traits::filetype, DB_READ )));
    else {
      // try open with DB_APPEND, if failed, do DBCreate
      // TODOL edge case: the file exists by accident so is garbage, but we need a new file
      // Temporarily disable error string to be printed
      auto* errfunc = DBErrfunc();
      int errlvl = DBErrlvl();
      DBShowErrors(DB_NONE, errfunc);
      DBfile* db = DBOpen( filename.c_str(), traits::filetype, DB_APPEND );
      DBShowErrors(errlvl, errfunc);
      if ( NULL == db ) {
        db = DBCreate( filename.c_str(), traits::create_mode, traits::create_target, NULL, traits::filetype );
      }
      dbfile.reset( new DBfileHandle(db));
    }
    return dbfile;
  }

}

namespace silo {
  std::string errmsg() {
    return std::string (DBErrFuncname()) + " : " + std::string (DBErrString());
  }
}
