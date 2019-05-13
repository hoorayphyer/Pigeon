#include "silopp/silo++.hpp"
#include <silo.h>
#include "mpipp/mpi++.hpp"
#include <pmpio.h> // needs mpi

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
  }

  template < Mode mode >
  file_t open(  std::string filename ) {
    file_t dbfile;

    if constexpr ( Mode::Read == mode )
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


  template file_t open<Mode::Read>( std::string filename );
  template file_t open<Mode::Write>( std::string filename );
}

namespace silo :: pmpio {
  void file_free( DBfileHandle* p ) {
    if ( p && (*p).baton_h ) {
      if ( (*p).file_h ) PMPIO_HandOffBaton( (*p).baton_h, (*p).file_h );
      PMPIO_Finish( (*p).baton_h );
    }
  }

  void* createCb( const char * fname, const char * dname, void * udata ) {
    DBfile * dbFile = DBCreate( fname, traits::create_mode, traits::create_target, NULL, traits::filetype );
    if ( dbFile ) {
      // create a directory in the silo file
      DBMkDir( dbFile, dname );
      // set current directory within the silo file
      DBSetDir( dbFile, dname );
    }

    return (void *) dbFile;
  }

  void* openCb( const char* fname, const char* nsname, PMPIO_iomode_t ioMode, void * udata ) {
    DBfile *siloFile = DBOpen(fname, traits::filetype, ioMode == PMPIO_WRITE ? DB_APPEND : DB_READ);
    if (siloFile)
      {
        if (ioMode == PMPIO_WRITE)
          DBMkDir(siloFile, nsname);
        DBSetDir(siloFile, nsname);
      }
    return (void *) siloFile;
  }

  template< Mode mode >
  file_t open( std::string filename, std::string dirname, const mpi::Comm& comm, int num_files ) {
    file_t dbfile;

    constexpr int mpi_tag = 147;

    auto pmpio_mode = Mode::Read == mode ? PMPIO_READ : PMPIO_WRITE;

    PMPIO_baton_t* baton_h = PMPIO_Init( num_files, pmpio_mode, comm, mpi_tag, createCb, openCb, PMPIO_DefaultClose, NULL );

    DBfile* file_h = (DBfile*) PMPIO_WaitForBaton( baton_h, filename.c_str(), dirname.c_str() );
    dbfile.reset( new DBfileHandle{ baton_h, file_h } );

    return dbfile;
  }

  template file_t open<Mode::Read>( std::string, std::string, const mpi::Comm&, int );
  template file_t open<Mode::Write>( std::string, std::string, const mpi::Comm&, int );

}

namespace silo {
  int to_group_rank( const mpi::Comm& comm, int num_files ) {
    // RATIONALE: rank / num_in_each_file, where num_in_each_file = ceiling( comm.size / num_files )
    return comm.rank() / ( comm.size() / num_files + ( comm.size() % num_files != 0 ) );
  }
}
