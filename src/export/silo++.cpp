#include "export/silo++.hpp"
#include "parallel/mpi++.hpp"

#include <mpi.h> // needed by pmpio
#include <silo.h>
#include <pmpio.h>
#include <memory>

namespace silo :: traits {
  bool display_guard = true;
  int create_mode = DB_CLOBBER;
  int create_target = DB_LOCAL;
  int filetype = DB_HDF5;
}

// NOTE the raw handle is DBfile* and PMPIO_baton_t*
// TODOL add is_detected check for raw_handle type, or is_same check

namespace silo {
  DBfile_t::~DBfile_t() {
    file_h.reset(); // file_h goes first
    _baton.reset();
  }

  // TODOL check correctness
  void close( DBfile_t dbfile ) {
    dbfile.file_h.reset();
    dbfile._baton.reset();
  }

  // TODOL check file existence
  template < Mode mode >
  DBfile_t open(  std::string filename ) {
    DBfile_t dbfile;

    DBfile* raw_file = nullptr;

    if constexpr ( Mode::Read == mode )
      raw_file = DBOpen( filename.c_str(), traits::filetype, DB_READ );
    else if ( Mode::Write == mode )
      // TODOL what if an existing file?
      raw_file = DBCreate( filename.c_str(), traits::create_mode, traits::create_target, NULL, traits::filetype );
    else
      raw_file = DBOpen( filename.c_str(), traits::filetype, DB_APPEND );

    dbfile.file_h = Handle( &raw_file, []( DBfile** f ) { if (f && *f) DBClose(*f); } );

    return dbfile;
  }


  template DBfile_t open<Mode::Read>( std::string filename );
  template DBfile_t open<Mode::Write>( std::string filename );
  template DBfile_t open<Mode::Append>( std::string filename );
}

namespace silo :: pmpio {
  int num_files = 2;

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
  DBfile_t open(  std::string dirname, const mpi::Comm& comm ) {
    static_assert( mode != Mode::Append, "PMPIO doesn't support Mode::Append");

    DBfile_t dbfile;

    constexpr int mpi_tag = 147;
    PMPIO_baton_t* raw_baton = PMPIO_Init( num_files, Mode::Read == mode ? PMPIO_READ : PMPIO_WRITE, comm.hdl, mpi_tag, createCb, openCb, PMPIO_DefaultClose, NULL );
    dbfile._baton = Handle( &raw_baton, []( PMPIO_baton_t** b ) { if (b && *b) PMPIO_Finish( *b ); } );

    auto filename = dirname + std::to_string( PMPIO_GroupRank( dbfile._baton, comm.rank() ) )+".silo";
    auto silo_dname = "proc" + std::to_string( comm.rank() );

    DBfile* raw_dbfile = (DBfile*) PMPIO_WaitForBaton( dbfile._baton, filename.c_str(), silo_dname.c_str() );
    dbfile.file_h = Handle( &raw_dbfile, [&baton = dbfile._baton]( DBfile** f) { if ( f && *f && baton != nullhdl ) PMPIO_HandOffBaton(baton, *f);} );

    return dbfile;
  }


  template DBfile_t open<Mode::Read>(  std::string dirname, const mpi::Comm& comm );
  template DBfile_t open<Mode::Write>(  std::string dirname, const mpi::Comm& comm );

}

// namespace silo{

//   std::function<void(DBfile*)> put_mesh;
//   std::function<void(DBfile*)> put_var;
//   std::function<void(DBfile*, int timestep)> put_multimesh;
//   std::function<void(DBfile*, int timestep, std::string name_sth, const int type)> put_multivar;
//   std::unique_ptr<PMPIO> pmpio;


//   template<CoordType Coord, int DIM>
//   auto set_putters( const Grid& grid, const int comm_size, const DBPane_DataExport& pane ) {

//     typename CoordToScales<Coord>::type scales;
//     constexpr bool linearMesh = (Coord == CoordType::CARTESIAN || Coord == CoordType::CYLINDRICAL);
//     // If a coord system is axisymmetric then we swap y and z axis
//     constexpr bool axisymmetry = decltype(scales)::axisymmetric;
//     constexpr char MESHNAME[] = "aperture_mesh";

//     std::array<int,DIM> dims;
//     for (int i = 0; i < DIM; ++i )
//       dims[i] = grid.dims[i];


//     put_mesh =
//       [=] (DBfile* dbfile, const Grid& grid, DBoptlist* optlist) {
//         constexpr auto COORDTYPE = (linearMesh ? DB_COLLINEAR : DB_NONCOLLINEAR);
//         // // FIXME TODO fix this and do we really need DB_NONCOLLINEAR??
//         // tempMesh.gridPoints.resize(dim);
//         // for (int i = 0; i < dim; i++) {
//         //   if (linearMesh)
//         //     tempMesh.gridPoints[i].resize(grid.dims[i], 0.0);
//         //   else
//         //     tempMesh.gridPoints[i].resize(grid.size(), 0.0);
//         //   tempMesh.gridArray[i] = tempMesh.gridPoints[i].data();
//         //   lowOffset[i] = grid.guard[i] - 1;
//         //   hiOffset[i] = grid.guard[i] - 1;
//         // }

//         // if (linearMesh) {
//         //   for (int i = 0; i < dim; i++) {
//         //     for (int j = 0; j < grid.dims[i]; ++j) {
//         //       tempMesh.gridPoints[i][j] = grid.pos(i, j, StaggerType::STAGGER_MINUS);
//         //     }
//         //   }
//         // } else {
//         //   for (int i = 0; i < dim; i++) {
//         //     for (int j = 0; j < grid.size(); ++j) {
//         //       Index idx(j, grid.extent());
//         //       Vec3<Scalar> pos(grid.pos(0, idx.x, StaggerType::STAGGER_MINUS),
//         //                        grid.pos(1, idx.y, StaggerType::STAGGER_MINUS),
//         //                        grid.pos(2, idx.z, StaggerType::STAGGER_MINUS));

//         //       scales.PosToCartesian(pos);
//         //       if (axisymmetry) {
//         //         std::swap(pos[1], pos[2]);
//         //       }

//         //       tempMesh.gridPoints[i][j] = pos[i];
//         //     }
//         //   }
//         // }


//         // DBPutQuadmesh(dbfile, MESHNAME, NULL, mesh.gridArray,
//         //               dims.data(), DIM, DB_FLOAT, COORDTYPE, optlist);
//       };

//     put_var =
//       [=] ( DBfile* dbfile, std::string name, const MultiArray<Scalar>& field ) {
//         //FIXME TODO do a downsampling here?
//         DBPutQuadvar1(dbfile, name.c_str(), MESHNAME, field.data(), dims.data(), DIM, NULL, 0, DB_FLOAT, DB_NODECENT, NULL); };
//   }

//   auto put_multi_sth =
//     [&pmpio=pmpio, nblocks=comm_size, isPmpio=pane.isPmpio, prefix=pane.prefix] (DBfile* dbfile, int timestep, std::string name_sth, const int type_sth, auto db_put_multisth) {
//       std::vector<char*> names(nblocks);

//       for (int i = 0; i < nblocks; i++) {
//         char* name = new char[100];
//         if ( isPmpio ) {
//           sprintf(name, "group%d/%s%06d.d:proc%d/%s", pmpio->GroupRank(i), prefix.c_str(), timestep, i, name_sth.c_str());
//         } else {
//           sprintf(name, "rank%d/%s%06d.d:%s", i, prefix.c_str(), timestep, name_sth.c_str());
//         }
//         names.push_back(name);
//       }

//       DBoptlist* optlist = DBMakeOptlist(1);
//       DBAddOption(optlist, DBOPT_MB_BLOCK_TYPE, &type_sth);
//       db_put_multisth(dbfile, name_sth.c_str(), nblocks, names.data(), NULL, optlist);
//       DBFreeOptlist(optlist);

//       for ( auto& elm : names )
//         delete[] elm;
//     };

//   put_multimesh = [=] (DBfile* dbfile, int timestep ) { return put_multi_sth(dbfile, timestep, MESHNAME, (linearMesh ? DB_QUAD_RECT : DB_QUAD_CURV), DBPutMultimesh); };

//   put_multivar = [=] (DBfile* dbfile, int timestep, std::string name_sth, const int type ) { return put_multi_sth(dbfile, timestep, name_sth, type, DBPutMultivar); };

// }
