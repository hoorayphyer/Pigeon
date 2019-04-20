#include "AperData_IO.h"
#include <silo.h>
#include "Types.h"
#include "Logger.h"
#include "InfoCollector.h"
#include "AperParams.h"
#include "mpi_wrapper.h"
#include "boost/filesystem.hpp"
#include <pmpio.h>

#include <type_traits>
#include <boost/fusion/adapted/struct/adapt_struct.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
// the following is also needed to compile boost::fusion::result_of::size<SinglePtc>, otherwise it complains about incomplete type
#include <boost/fusion/sequence.hpp>

#include "ProcessRoster.h"

// ----- MetaFunctions for Write/Read particles -----

template < typename T >
inline int ToDBTypeImpl();

template <>
inline int ToDBTypeImpl<int>() { return DB_INT; }

template <>
inline int ToDBTypeImpl<unsigned int>() { return DB_INT; }

template <>
inline int ToDBTypeImpl<double>() { return DB_DOUBLE; }

template <>
inline int ToDBTypeImpl<float>() { return DB_FLOAT; }

template <>
inline int ToDBTypeImpl<char>() { return DB_CHAR; }

// a wrapper over ToDBTypeImpl that strips off all const qualifiers and references
template < typename T >
inline int ToDBType() {
  return ToDBTypeImpl< typename std::remove_const< typename std::remove_reference<T>::type >::type >();
}

BOOST_FUSION_ADAPT_STRUCT( QCarrier,
                           x, dx, p, cell, flag, Rc, track_id )

BOOST_FUSION_ADAPT_STRUCT( QNeutral,
                           x, p, cell, flag, path_left, E, track_id )

// define some type traits for the elements of a SinglePtc
template <typename T>
struct TypeTraits {
  using base_type = T;
  static const int num_base_types = 1;
  static const base_type& get( const T& a, int idx ) { return a; }
  static base_type& get( T& a, int idx ) { return a; }
};

// specialization for the Vec3 class
template <typename T>
struct TypeTraits< Vec3<T> > {
  using base_type = T;
  static const int num_base_types = VECTOR_DIM;
  static const base_type& get( const Vec3<T>& a, int idx ) { return a[idx]; }
  static base_type& get( Vec3<T>& a, int idx ) { return a[idx]; }
};

// specialization for std::array
template <typename T, std::size_t N>
struct TypeTraits< std::array<T,N> > {
  using base_type = T;
  static const int num_base_types = N;
  static const base_type& get( const std::array<T,N>& a, int idx ) { return a[idx]; }
  static base_type& get( std::array<T,N>& a, int idx ) { return a[idx]; }
};

template < typename SinglePtc, int I_Elem, bool Cond >
struct PtcLoop {
  // NOTE at_c returns a type reference
  using ElemT = typename std::remove_reference< typename boost::fusion::result_of::at_c<SinglePtc, I_Elem>::type >::type;
  using BaseElemT = typename TypeTraits<ElemT>::base_type;

  // the return type is const char *
  static constexpr decltype(auto) ElementName() {
    return boost::fusion::extension::struct_member_name<SinglePtc, I_Elem>::call();
  }

  template < typename Filter >
  static void Write ( DBfile* dbfile, const Particles<SinglePtc>& particles, void* memory, int num, const Filter& filter ) {
    constexpr int NDIMS = 1;
    constexpr int num_base_types = TypeTraits<ElemT>::num_base_types;
    int dims[NDIMS] = { num * num_base_types };

    auto ptr = reinterpret_cast< BaseElemT* >( memory );
    int idx_selected = 0;
    // NOTE the loop is over particles.Number() as opposed to num due to the use of filter
    for ( int i = 0; i < particles.Number(); ++i ) {
      if ( particles.IsEmpty(i) || !filter(particles, i) ) continue;
      const auto& elm_value = boost::fusion::at_c< I_Elem >( particles.PtcData()[i] );
      for ( int j = 0; j < num_base_types; ++j )
        ptr[ idx_selected * num_base_types + j ] = TypeTraits<ElemT>::get( elm_value, j );
      ++idx_selected;
    }

    DBWrite( dbfile, ElementName(), ptr, dims, NDIMS, ToDBType<BaseElemT>() );

    // iterate
    PtcLoop< SinglePtc, I_Elem + 1, ( I_Elem + 1 ) < boost::fusion::result_of::size<SinglePtc>::value >::Write( dbfile, particles, memory, num, filter );
  }

  static void Read ( DBfile* dbfile, Particles<SinglePtc>& particles, void* memory, const int& num, const int& offset_from, const int& offset_to ) {
    constexpr int NDIMS = 1;
    constexpr int num_base_types = TypeTraits<ElemT>::num_base_types;
    // define some helper variables for DBReadVarSlice
    int db_offset[NDIMS] = { offset_from * num_base_types };
    int db_length[NDIMS] = { num * num_base_types };
    int db_stride[NDIMS] = {1};
    int db_ndims = 1;

    auto ptr = reinterpret_cast< BaseElemT* >( memory );
    DBReadVarSlice( dbfile, ElementName(), db_offset, db_length, db_stride, db_ndims, ptr );
    for ( int i = 0; i < num; ++i ) {
      auto& elm_value = boost::fusion::at_c< I_Elem >( particles.PtcData()[ offset_to + i ] );
      for ( int j = 0; j < num_base_types; ++j )
        TypeTraits<ElemT>::get(elm_value, j) = ptr[ i * num_base_types + j ];
    }

    // iterate
    PtcLoop< SinglePtc, I_Elem + 1, ( I_Elem + 1 ) < boost::fusion::result_of::size<SinglePtc>::value >::Read( dbfile, particles, memory, num, offset_from, offset_to );
  }

};

// specialize for Cond = false
template < typename SinglePtc, int I >
struct PtcLoop < SinglePtc, I, false > {
  template <typename Filter>
  static void Write ( DBfile* dbfile, const Particles<SinglePtc>& particles, void* memory, int num, const Filter& filter ) {}

  // do some wrap up work
  static void Read ( DBfile* dbfile, Particles<SinglePtc>& particles, void* memory, const int& num, const int& offset_from, const int& offset_to ) {
    particles.SetNum( particles.Number() + num );
  }

};

// ----- Pmpio Callbacks Adapted From <pmpio> -----

void * PMPIO_CustomCreate(const char *fname, const char *nsname, void *userData)
{
  DBfile *siloFile = DBCreate(fname, DB_CLOBBER, DB_LOCAL, "PMPIO_CustomCreate", DB_HDF5);
  if (siloFile) {
    if ( !DBInqVarExists( siloFile, nsname ) ) {
      DBMkDir(siloFile, nsname);
    }
    DBSetDir(siloFile, nsname);
  }
  return (void *) siloFile;
}

void * PMPIO_CustomOpen(const char *fname, const char *nsname, PMPIO_iomode_t ioMode, void *userData)
{
  DBfile *siloFile = DBOpen(fname, DB_HDF5,
                            ioMode == PMPIO_WRITE ? DB_APPEND : DB_READ);
  if (siloFile) {
    if ( ioMode == PMPIO_WRITE && !DBInqVarExists( siloFile, nsname ) )
      DBMkDir(siloFile, nsname);
    DBSetDir(siloFile, nsname);
  }
  return (void *) siloFile;
}



// ----- Implementation of Class Member Functions -----

void AperData_IO::WriteSnapshot( Scalar time, int timestep, const AperParams& params, const AperData& data, const MPICommunicator& comm, const ProcessRoster* roster, std::string snapshot_dir, const int num_pmpio_files ) const {
  //------ generate the communicator of active processes ( primary or replica ) for pmpio
  // NOTE: one could also just use MPI_COMM_WORLD here to init pmpio and then let only active processes output, which could avoid generating the active_comm. The downside of this is that the inclusion of idles may break the even distribution of active processes in each pmpio silo file.
  auto t0 = high_resolution_clock::now();
  std::vector<int> active_ranks;
  const auto& world = comm.world();
  if ( world.is_root() ) {
    active_ranks = roster->active_ranks();
    int size = active_ranks.size();
    world.broadcast( &size, 1, world.root() );
    world.broadcast( active_ranks.data(), size, world.root() );
  } else {
    int size = 0;
    world.broadcast( &size, 1, world.root() );
    active_ranks.resize( size );
    world.broadcast( active_ranks.data(), size, world.root() );
  }

  // idles can now exit
  if ( comm.is_idle() ) return;
  auto t1 = high_resolution_clock::now();
  auto dur = duration_cast<clock_cycle>(t1 - t0);
  InfoCollector::Instance().contents.t_saveSnapshot_activeGroup = dur.count();

  MPI_Group grp_active;
  MPI_Group_incl( world.group(), active_ranks.size(), active_ranks.data(), &grp_active );
  MPI_Comm comm_active;
  MPI_Comm_create_group(world.comm(), grp_active, 147, &comm_active);
  MPI_Group_free( &grp_active );
  int comm_active_rank = 0;
  MPI_Comm_rank(comm_active, &comm_active_rank);

  //------ prepare pmpio ------
  const auto& ensemble = comm.ensemble();

  PMPIO_baton_t* pmpio_bat = PMPIO_Init( num_pmpio_files, PMPIO_WRITE, comm_active, 147, PMPIO_CustomCreate, PMPIO_CustomOpen, PMPIO_DefaultClose, NULL );

  // assume snapshot_dir ends with "/"
  std::string snapshot_filename = snapshot_dir + "snapshot_timestep" + std::to_string(timestep) + "_group" + std::to_string( PMPIO_GroupRank(pmpio_bat, comm_active_rank ) ) + ".silo";

  std::string silo_dirname( "/ensemble_" + std::to_string(ensemble.label() ) + "/" );

  auto dbfile = (DBfile*) PMPIO_WaitForBaton( pmpio_bat, snapshot_filename.c_str(), silo_dirname.c_str() );

  //------ Write data ------
  // write global variables
  if ( !DBInqVarExists( dbfile, "/timestep") ) {
    int dims[1] = {1};
    DBWrite( dbfile, "/timestep", &timestep, dims, 1, ToDBType<decltype(timestep)>() );
  }
  if ( !DBInqVarExists( dbfile, "/num_ensembles") && !comm.cartesian().isNull() ) {
    int dims[1] = {1};
    int num_ensembles = comm.cartesian().size();
    DBWrite( dbfile, "/num_ensembles", &num_ensembles, dims, 1, ToDBType<decltype(num_ensembles)>() );
  }

  // write ensemble-common variables
  if ( !DBInqVarExists( dbfile, (silo_dirname + "ens_size").c_str() ) ) {
    { // write ens size
      int dims[1] = {1};
      int ens_size = ensemble.size();
      DBWrite( dbfile, (silo_dirname + "ens_size").c_str(), &ens_size, dims, 1, ToDBType<decltype(ens_size)>());
    }
    { // write ens label. Explicitly save this to avoid the use of regex to extract ens label from the directory name
      int dims[1] = {1};
      int ens_label = ensemble.label();
      DBWrite( dbfile, (silo_dirname + "ens_label").c_str(), &ens_label, dims, 1, ToDBType<decltype(ens_label)>());
    }
    { // write cartesian position
      int cart_pos[3] = { params.ens_specs.coordinate[0], params.ens_specs.coordinate[1], params.ens_specs.coordinate[2] };
      int dims[1] = {3};
      DBWrite( dbfile, (silo_dirname + "cart_pos").c_str(), cart_pos, dims, 1, ToDBType<decltype(cart_pos[0])>());
    }
  }

  { // write the tracker counter statuses under the directory /trackers. The
    // name is the world rank, and the content is an array of current counter
    // values in the order of electron, positron, ion and photon.
    if ( !DBInqVarExists( dbfile, "/trackers/" ) ) {
      DBMkDir( dbfile, "/trackers/" );
    }
    std::array< unsigned int, NUM_PTC_SPECIES + 1 > counters;
    for ( int i = 0; i < NUM_PTC_SPECIES; ++i ) {
      auto ptcType = static_cast<ParticleType>(i);
      if ( data.particles.find(ptcType) != data.particles.end() ) {
        counters[i] = data.particles.at(ptcType).Tracker().ReadCounter();
      } else {
        counters[i] = 0;
      }
    }
    counters[static_cast<int>(ParticleType::PHOTON)] = data.photons.Tracker().ReadCounter();
    int dims[1] = { counters.size() };
    DBWrite(dbfile, ("/trackers/world_rank" + std::to_string(world.rank()) ).c_str(), counters.data(), dims, 1, ToDBType<decltype(counters[0])>() );
  }

  // world root writes the current ensemble size map for the sake of restarts
  if ( world.is_root() ) {
    const auto& label2ranks = roster->label2ranks;
    std::vector<int> size_map( label2ranks.size() );
    for ( unsigned int i = 0; i < size_map.size(); ++i )
      size_map[i] = label2ranks[i].size();

    int dims[1] = { size_map.size() };
    DBWrite( dbfile, "/size_map", size_map.data(), dims, 1, ToDBType<decltype(size_map[0])>() );
  }

  // make separate directory for primaries and replicas
  std::string work_dirname (  comm.is_primary() ? "primary/" :  "replica" + std::to_string(ensemble.rank()) + "/" );
  DBMkDir( dbfile, work_dirname.c_str() );
  DBSetDir( dbfile, work_dirname.c_str() );

  for ( const auto& elm : data.particles ) {
    const auto& particles = elm.second;
    int num = particles.Number();
    if ( 0 == num ) continue;
    WriteParticle( dbfile, particles, []( auto&&... args ) { return true; } );
  }
  if ( data.photons.Number() > 0 )
    WriteParticle( dbfile, data.photons, []( auto&&... args ) { return true; } );

  if ( comm.is_primary() ) {
    WriteField( dbfile, data.Efield, "E");
    WriteField( dbfile, data.Bfield, "B");
    WriteField( dbfile, data.Jfield, "J");
    WriteField( dbfile, data.pairCreationEvents, "PairCreation");

    WriteExtra( dbfile, InfoCollector::Instance().ssProxy );
  }

  //------ Wrap up ------
  PMPIO_HandOffBaton( pmpio_bat, dbfile );
  PMPIO_Finish( pmpio_bat );
  MPI_Comm_free( &comm_active );
}

// NOTE it is assumed that AperData is initialized after calling Init
// snapshot_dir is the directory that holds all silo files. A typical value is Data(date)/snapshots/snapshot_timestep(number)
void AperData_IO::ReadSnapshot(  AperData& data, const AperParams& params, const MPICommunicator& comm, std::string snapshot_dir ) const {
  if ( !boost::filesystem::exists( snapshot_dir ) ) {
    throw std::runtime_error("ERROR in ReadSnapshot : " + snapshot_dir + " doesn't exist!");
  }

  boost::filesystem::path snapshot_dir_path(snapshot_dir);
  using dir_itr_type = boost::filesystem::directory_iterator;

  // assume the cartesian topology is the same between restarts
  const auto& pos = params.ens_specs.coordinate;

  // erase particle arrays first
  for ( auto& elm : data.particles ) {
    auto& particles = elm.second;
    particles.Erase( 0, particles.NumMax() );
  }
  data.photons.Erase( 0, data.photons.NumMax() );

  if ( comm.is_idle() ) return;
  const auto& ensemble = comm.ensemble();

  // loop over all snapshot .silo files
  ParticleArrayPartitioner partitioner( ensemble.size() );
  for ( dir_itr_type it(snapshot_dir_path); it != dir_itr_type(); ++it ) {
    DBfile* dbfile = DBOpen( it->path().string().c_str(), DB_HDF5, DB_READ );

    { // read tracker if applicable
      int world_rank = comm.world().rank();
      std::string tracker_var = "/trackers/world_rank" + std::to_string(world_rank);
      if ( DBInqVarExists( dbfile, tracker_var.c_str() ) ) {
        std::array<unsigned int, NUM_PTC_SPECIES + 1> counters;
        DBReadVar( dbfile, tracker_var.c_str(), counters.data() );
        for ( int i = 0; i < NUM_PTC_SPECIES; ++i ) {
          auto ptcType = static_cast<ParticleType>(i);
          if ( data.particles.find(ptcType) != data.particles.end() ) {
            data.particles.at(ptcType).Tracker().SetTrackCounter( counters[i] );
          }
        }
        // photon
        data.photons.Tracker().SetTrackCounter( counters[static_cast<int>(ParticleType::PHOTON)] );
      }
    }

    // find the only relevant ensemble/ directory if exists
    std::string my_ensemble_dir = "/ensemble_" + std::to_string( ensemble.label() ) + "/";
    if ( DBInqVarExists( dbfile, my_ensemble_dir.c_str() ) ) {

      DBSetDir( dbfile, my_ensemble_dir.c_str() );
      DBtoc* dbtoc_ens = DBGetToc( dbfile );
      // loop over primary/ and replica/ of current pos/
      for ( int j = 0; j < dbtoc_ens->ndir; ++j ) {
        // switch to primary/ or replica/
        DBSetDir( dbfile, dbtoc_ens->dir_names[j] );

        // read particles
        for ( auto& elm : data.particles ) {
          ReadParticle( dbfile, elm.second, ensemble.rank(), partitioner );
        }

        // read photons
        ReadParticle( dbfile, data.photons, ensemble.rank(), partitioner );

        // switch back
        DBSetDir( dbfile, ".." );
        // restore toc because it might be destroyed by any DBSetDir in previous lines
        dbtoc_ens = DBGetToc( dbfile );

        partitioner.AdvanceShift();
      }

      // primary reads from ensemble/primary for fields if applicable
      if ( DBInqVarExists(dbfile, "primary/") && comm.is_primary() ) {
        DBSetDir( dbfile, "primary");

        ReadField( dbfile, data.Efield, "E" );
        ReadField( dbfile, data.Bfield, "B" );
        ReadField( dbfile, data.Jfield, "J" );
        ReadField( dbfile, data.pairCreationEvents, "PairCreation");

        ReadExtra( dbfile, InfoCollector::Instance().ssProxy );

        DBSetDir( dbfile, ".." );
      }

      DBSetDir( dbfile, ".." );
    }

    DBClose( dbfile );
  }

}

void AperData_IO::WriteTrackingData( int timestep, AperData& data, const MPIWorldCommunicator& world, const Grid& grid, std::string track_dir, const int num_pmpio_files ) const {
  // use separate silo files for each timestep in order to be restart-friendly
  // each process will write tracked particles data and tracked pair creation events if any.
  // in-file structure is : /pair_creation_tracking/data_named_with_world_rank, /tracked_particles/species/particle_data_named_with_world_rank
  PMPIO_baton_t* pmpio_bat = PMPIO_Init( num_pmpio_files, PMPIO_WRITE, world.comm(), 147, PMPIO_CustomCreate, PMPIO_CustomOpen, PMPIO_DefaultClose, NULL );

  // assume track_dir ends with "/"
  std::string pc_track_filename = track_dir + "tracking_timestep" + std::to_string(timestep) + "_group" + std::to_string( PMPIO_GroupRank( pmpio_bat, world.rank() ) ) + ".silo";

  std::string silo_dirname( "/" );
  auto dbfile = (DBfile*) PMPIO_WaitForBaton( pmpio_bat, pc_track_filename.c_str(), silo_dirname.c_str() );

  // ----- save tracked particles data -----
  WriteTrackedParticles( dbfile, data, grid, track_dir, world.rank() );

  // ----- save pair creation data -----
  WriteTrackedPairCreation( dbfile, data.pairCreationTracker, track_dir, world.rank() );

  //------ Wrap up ------
  PMPIO_HandOffBaton( pmpio_bat, dbfile );
  PMPIO_Finish( pmpio_bat );
}

int AperData_IO::ReadTimestep( std::string snapshot_dir ) const {
  if ( !boost::filesystem::exists( snapshot_dir ) ) {
    throw std::runtime_error("ERROR in ReadTime : " + snapshot_dir + " doesn't exist!");
  }

  boost::filesystem::path snapshot_dir_path(snapshot_dir);
  boost::filesystem::directory_iterator dir_it(snapshot_dir_path);

  // NOTE: dir_it points to the begining file in the directory
  int timestep = 0.0;
  DBfile* dbfile = DBOpen(dir_it->path().string().c_str(), DB_HDF5, DB_READ);
  DBReadVar(dbfile, "timestep", &timestep);
  DBClose( dbfile );

  return timestep;
}

std::vector<int> AperData_IO::ReadSizeMap( std::string snapshot_dir ) const {
  std::vector<int> size_map;

  if ( !boost::filesystem::exists( snapshot_dir ) ) {
    throw std::runtime_error("ERROR in ReadSizeMap : " + snapshot_dir + " doesn't exist!");
  }

  boost::filesystem::path snapshot_dir_path(snapshot_dir);
  using dir_itr_type = boost::filesystem::directory_iterator;

  // loop over all snapshot .silo files
  for ( dir_itr_type it(snapshot_dir_path); it != dir_itr_type(); ++it ) {
    DBfile* dbfile = DBOpen(it->path().string().c_str(), DB_HDF5, DB_READ);
    if ( DBInqVarExists(dbfile, "/size_map") ) {
      size_map.resize( DBGetVarLength( dbfile, "/size_map" ) );
      DBReadVar( dbfile, "/size_map", size_map.data() );
      DBClose( dbfile );
      break;
    }
    DBClose( dbfile );
  }

  return size_map;
}

std::vector<size_t> AperData_IO::ReadParticleLoadMap( std::string snapshot_dir, size_t (*load_calculator)(size_t N_e, size_t N_p, size_t N_i, size_t N_ph) ) const {
  if ( !boost::filesystem::exists( snapshot_dir ) ) {
    throw std::runtime_error("ERROR in ReadParticleLoadMap : " + snapshot_dir + " doesn't exist!");
  }

  boost::filesystem::path snapshot_dir_path(snapshot_dir);
  using dir_itr_type = boost::filesystem::directory_iterator;

  int num_ensembles = 0;
  // read num_ensembles
  for ( dir_itr_type it(snapshot_dir_path); it != dir_itr_type(); ++it ) {
    DBfile* dbfile = DBOpen(it->path().string().c_str(), DB_HDF5, DB_READ);
    if ( DBInqVarExists(dbfile, "/num_ensembles") ) {
      DBReadVar( dbfile, "/num_ensembles", &num_ensembles );
      DBClose( dbfile );
      break;
    }
    DBClose( dbfile );
  }

  if ( 0 == num_ensembles ) {
    throw std::runtime_error("Fail to read num_ensembles in ReadParticleLoadMap");
  }

  std::vector<size_t> ptc_load_map( num_ensembles );
  std::fill( ptc_load_map.begin(), ptc_load_map.end(), 0 );

  // loop over all snapshot .silo files
  for ( dir_itr_type it(snapshot_dir_path); it != dir_itr_type(); ++it ) {
    DBfile* dbfile = DBOpen(it->path().string().c_str(), DB_HDF5, DB_READ);
    ReadParticleLoadMapFromOneFile( ptc_load_map, dbfile, load_calculator );
    DBClose( dbfile );
  }

  return ptc_load_map;
}

// The Filter function should support operator( const Particles<SinglePtc>& , int idx ).
template < typename SinglePtc, typename Filter >
void AperData_IO::WriteParticle( DBfile* dbfile, const Particles<SinglePtc>& particles, Filter filter ) const {
  int num = 0;
  for ( int i = 0; i < particles.Number(); ++i ) {
    if ( !particles.IsEmpty(i) && filter(particles, i) ) ++num;
  }
  // if there is no particle, skip write
  if ( 0 == num ) {
    Logger::print(Logger::gVerbosityLvl, "Zero", particles.NameStr(), "selected in array. Skip Writing");
    return;
  }

  std::string species_dir = particles.NameStr();
  // mkdir with name of particle species
  DBMkDir( dbfile, species_dir.c_str() );
  DBSetDir( dbfile, species_dir.c_str() );

  std::unique_ptr<double[]> memory ( new double [num * VECTOR_DIM] );
  // WriteParticleImpl( dbfile, particles, memory.get(), num );
  PtcLoop<SinglePtc, 0, 0 < boost::fusion::result_of::size<SinglePtc>::value >::Write( dbfile, particles, memory.get(), num, filter );

  DBSetDir( dbfile, ".." );
}

template <typename SinglePtc, typename Partitioner >
void AperData_IO::ReadParticle( DBfile* dbfile, Particles<SinglePtc>& particles, int ens_rank, const Partitioner& partitioner ) const {
  std::string species_dir = particles.NameStr();
  // check if the directory of the species exists. If not, return
  if ( !DBInqVarExists( dbfile, species_dir.c_str() ) ) {
    Logger::print_screen(Logger::gVerbosityLvl, "Rank", Logger::thisRank, ",", species_dir, "doesn't exist");
    return;
  }

  int num_tot = DBGetVarLength( dbfile, (species_dir + "/cell").c_str() );
  int offset_from = 0;
  int num = 0;
  partitioner( num_tot, ens_rank, offset_from, num );
  int offset_to = particles.Number();
  // Logger::print_screen(Logger::gVerbosityLvl, "Rank", Logger::thisRank, "reads in", num, particles.NameStr() );

  // allocate a large enough memory
  std::unique_ptr<double> memory( new double [ num * VECTOR_DIM ] );

  DBSetDir( dbfile, species_dir.c_str() );
  // ReadParticleImpl( dbfile, particles, (void*) memory.get(), num, offset_from, offset_to );
  PtcLoop<SinglePtc, 0, (0 < boost::fusion::result_of::size<SinglePtc>::value) >::Read( dbfile, particles, (void*) memory.get(), num, offset_from, offset_to );
  DBSetDir( dbfile, "..");

}

template < typename T, template <typename> class Field >
struct FieldTraits;

template < typename T >
struct FieldTraits< T, ScalarField > {
  static const int num_scalar_fields = 1;
  static T* get_ptr( ScalarField<T>& field, int idx ) { return field.ptr(); }
  static const T* get_ptr( const ScalarField<T>& field, int idx ) { return field.ptr(); }
};

template < typename T >
struct FieldTraits< T, VectorField > {
  static const int num_scalar_fields = 3;
  static T* get_ptr( VectorField<T>& field, int idx ) { return field.ptr(idx); }
  static const T* get_ptr( const VectorField<T>& field, int idx ) { return field.ptr(idx); }
};


template <typename T, template <typename> class Field >
void AperData_IO::WriteField( DBfile* dbfile, const Field<T>& field, std::string fieldname ) const {
  using Traits = FieldTraits<T, Field >;
  int db_datatype = ToDBType<T>();
  constexpr int ndims = 3;
  const auto& grid = field.grid();
  int dims[ndims] = { grid.dims[0], grid.dims[1], grid.dims[2] };

  for ( int i = 0; i < Traits::num_scalar_fields; ++i ) {
    std::string varname = (Traits::num_scalar_fields == 1) ? fieldname : fieldname + std::to_string(i+1);
    DBWrite( dbfile, varname.c_str(), Traits::get_ptr(field, i), dims, ndims, db_datatype );
  }

}

template <typename T, template <typename> class Field >
void AperData_IO::ReadField( DBfile* dbfile, Field<T>& field, std::string fieldname ) const {
  using Traits = FieldTraits<T, Field>;

  for ( int i = 0; i < Traits::num_scalar_fields; ++i ) {
    std::string varname = (Traits::num_scalar_fields == 1) ? fieldname : fieldname + std::to_string(i+1);
    DBReadVar( dbfile, varname.c_str(), Traits::get_ptr(field, i) );
  }
}

template < typename SinglePtc >
void WriteTrackedParticlesImpl( DBfile* dbfile, const Particles<SinglePtc>& particles, const Grid& grid, const std::string& save_dir ) {
  // will write particle's physical position, track_id, momentum
  int num_tracked = 0;
  for ( int i = 0; i < particles.Number(); ++i ) {
    if ( particles.IsEmpty(i) ) continue;
    const auto& ptc = particles.PtcData()[i];
    // don't save particles outside r = 12
    Scalar ln_r = grid.pos( 0, grid.getC1(ptc.cell), 0, ptc.x.x );
    if ( ln_r > 2.48 ) continue;
    if ( ParticleTracker::IsTracked( ptc.track_id ) ) ++num_tracked;
  }
  if ( 0 == num_tracked ) return;

  bool is_secondary = ParticleType::ELECTRON == particles.Attributes().ptcType || ParticleType::POSITRON == particles.Attributes().ptcType;

  std::unique_ptr<int[]> id_array ( new int [ 2 * num_tracked ] );
  std::unique_ptr<Scalar[]> pos_array ( new Scalar [ 3 * num_tracked ] );
  std::unique_ptr<Scalar[]> mom_array ( new Scalar [ 3 * num_tracked ] );

  std::unique_ptr<int[]> sec_array;
  if ( is_secondary ) {
    sec_array.reset( new int [num_tracked] );
  }

  int track_idx = 0;
  for ( int i = 0; i < particles.Number(); ++i ) {
    if ( particles.IsEmpty(i) ) continue;
    const auto& ptc = particles.PtcData()[i];
    if ( !ParticleTracker::IsTracked( ptc.track_id ) ) continue;
    // don't save particles outside r = 12
    Scalar ln_r = grid.pos( 0, grid.getC1(ptc.cell), 0, ptc.x.x );
    if ( ln_r > 2.48 ) continue;

    for ( int j = 0; j < 2; ++j ) {
      id_array[ track_idx * 2 + j ] = ptc.track_id[j];
    }
    auto pos = grid.pos_particle( ptc.cell, ptc.x );
    for ( int j = 0; j < 3; ++j ) {
      pos_array[ track_idx * 3 + j ] = pos[j];
      mom_array[ track_idx * 3 + j ] = ptc.p[j];
    }
    if ( is_secondary ) {
      sec_array[ track_idx ] = static_cast<int>( check_bit( ptc.flag, ParticleFlag::secondary ) );
    }
    ++track_idx;
  }

  if ( !DBInqVarExists( dbfile, save_dir.c_str()) )
    DBMkDir( dbfile, save_dir.c_str() );
  DBSetDir( dbfile, save_dir.c_str() );

  {
    int dims[1] = { 2 * num_tracked };
    DBWrite( dbfile, "track_id", id_array.get(), dims, 1, ToDBType<decltype(id_array[0])>() );
  }
  {
    int dims[1] = { 3 * num_tracked };
    DBWrite( dbfile, "position", pos_array.get(), dims, 1, ToDBType<decltype(pos_array[0])>() );
    DBWrite( dbfile, "momentum", mom_array.get(), dims, 1, ToDBType<decltype(mom_array[0])>() );
  }
  if ( is_secondary ) {
    int dims[1] = { num_tracked };
    DBWrite( dbfile, "is_secondary", sec_array.get(), dims, 1, ToDBType<decltype(sec_array[0])>() );
  }

  DBSetDir( dbfile, ".." );
}

void AperData_IO::WriteTrackedParticles( DBfile* dbfile, const AperData& data, const Grid& grid, const std::string& track_dir, int world_rank ) const {
  std::string ptc_tracking_dir = "tracked_particles/";
  if ( !DBInqVarExists( dbfile, ptc_tracking_dir.c_str() ) ) {
    DBMkDir( dbfile, ptc_tracking_dir.c_str() );
  }
  DBSetDir( dbfile, ptc_tracking_dir.c_str() );

  // deal with particles and photons together
  std::vector<ParticleType> ptcTypes;
  for ( auto it = data.particles.cbegin(); it != data.particles.cend(); ++it ) {
    ptcTypes.push_back(it->first);
  }
  ptcTypes.push_back(ParticleType::PHOTON);
  for ( auto ptcType : ptcTypes ) {
    std::string species_dir = PtcType2Str(ptcType) + "/";
    if ( !DBInqVarExists( dbfile, species_dir.c_str() ) ) {
      DBMkDir( dbfile, species_dir.c_str() );
    }
    DBSetDir( dbfile, species_dir.c_str() );
    if ( ptcType != ParticleType::PHOTON ) {
      // ptcType is guaranteed to exist, due to the way ptcTypes is created
      const auto& particles = data.particles.at(ptcType);
      WriteTrackedParticlesImpl( dbfile, particles, grid, "world_rank" + std::to_string(world_rank) );
    } else {
      const auto& photons = data.photons;
      WriteTrackedParticlesImpl( dbfile, photons, grid, "world_rank" + std::to_string(world_rank) );
    }
    DBSetDir( dbfile, ".." );
  }

  DBSetDir( dbfile, ".." );
}

void AperData_IO::WriteTrackedPairCreation( DBfile* dbfile, PairCreationTracker& tracker, const std::string& track_dir, int world_rank  ) const {
  std::string pc_tracking_dir = "pair_creation_tracking/";
  if ( !DBInqVarExists( dbfile, pc_tracking_dir.c_str() ) ) {
    DBMkDir( dbfile, pc_tracking_dir.c_str() );
  }
  DBSetDir( dbfile, pc_tracking_dir.c_str() );
  // only processes with non empty pc tracking data will save
  const auto& track_data = tracker.GetData();
  if ( track_data.size() != 0 ) {
    std::string varname = "world_rank" + std::to_string( world_rank );
    int dims[1] = { track_data.size() };
    DBWrite( dbfile, varname.c_str(), track_data.data(), dims, 1, ToDBType<decltype(track_data[0])>() );
  }
  // clear out tracker data
  tracker.ClearData();
}

void AperData_IO::WriteExtra( DBfile* dbfile, const SaveSnapshotProxy& proxy ) const {
  // // Write injection residuals
  // for( int b = 0; b < NUM_BOUNDARIES; ++b ) {
  //   const auto& inj_res = proxy.injection_residuals[b];
  //   if ( inj_res.size() == 0 ) continue;

  //   std::string varname("injection_residual" + std::to_string(b) );
  //   constexpr int NDIMS = 1;
  //   int dims[NDIMS] = { inj_res.size() };
  //   DBWrite( dbfile, varname.c_str(), inj_res.data(), dims, NDIMS, ToDBType<decltype(*inj_res.data())>() );
  // }

  // // Write injection residuals rng state
  // if ( 0 != proxy.rng_state_inj.size() ) {
  //   int ndims = 1;
  //   int dims[1] = { proxy.rng_state_inj.size() };
  //   DBWrite( dbfile, "rng_state_inj", proxy.rng_state_inj.data(), dims, ndims, ToDBType<decltype(*proxy.rng_state_inj.data())>() );
  // }


  // Write damping background
  if ( proxy.fBC_damping_E_bg.size() != 0 ) {

    if ( proxy.fBC_damping_E_bg.size() != proxy.fBC_damping_B_bg.size() ) {
      throw std::runtime_error("Sizes differ between E and B damping background in SSProxy when calling WriteExtra");
    }

    int ndims = 1;
    int dims[1] = { proxy.fBC_damping_E_bg.size() };
    DBWrite( dbfile, "damping_E_bg", proxy.fBC_damping_E_bg.data(), dims, ndims, ToDBType<decltype(*proxy.fBC_damping_E_bg.data())>() );
    DBWrite( dbfile, "damping_B_bg", proxy.fBC_damping_B_bg.data(), dims, ndims, ToDBType<decltype(*proxy.fBC_damping_B_bg.data())>() );

  }

  // // Write rad transfer rng state
  // if ( 0 != proxy.rng_state_rad.size() ) {
  //   int ndims = 1;
  //   int dims_rad[1] = { proxy.rng_state_rad.size() };
  //   DBWrite( dbfile, "rng_state_rad", proxy.rng_state_rad.data(), dims_rad, ndims, ToDBType<decltype(*proxy.rng_state_rad.data())>() );
  // }

}

// In reading rng states, std::string does not have non-const data() member. Hence this function.
// void ReadToString( DBfile* dbfile, const char * varname, std::string& dest ) {
//   int size = 0;

//   if ( DBInqVarExists( dbfile, varname ) ) {
//     size = DBGetVarLength( dbfile, varname );
//   }

//   if ( 0 == size ) {
//     return;
//   }

//   char * buf = new char[size];
//   DBReadVar( dbfile, varname, buf );

//   dest = std::string(buf, size);
//   delete [] buf;
// }

void AperData_IO::ReadExtra( DBfile* dbfile, SaveSnapshotProxy& proxy ) const {
  // // Read injection residuals and rng state
  // for ( int b = 0; b < NUM_BOUNDARIES; ++b ) {
  //   std::string varname("injection_residual" + std::to_string(b) );
  //   if ( DBInqVarExists( dbfile, varname.c_str()) == 0 ) continue;

  //   int size = DBGetVarLength( dbfile, varname.c_str() );
  //   proxy.injection_residuals[b].resize( size );
  //   DBReadVar(dbfile, varname.c_str(), proxy.injection_residuals[b].data());

  // }
  // // ReadToString( dbfile, "rng_state_inj", proxy.rng_state_inj );
  // proxy.isReadInjection = true;

  // Read damping background
  if ( DBInqVarExists( dbfile, "damping_E_bg") ) {
      int size = DBGetVarLength( dbfile, "damping_E_bg" );

      proxy.fBC_damping_E_bg.resize( size );
      proxy.fBC_damping_B_bg.resize( size );
      DBReadVar(dbfile, "damping_E_bg", proxy.fBC_damping_E_bg.data());
      DBReadVar(dbfile, "damping_B_bg", proxy.fBC_damping_B_bg.data());

      proxy.isReadDampingBg = true;
  }

  // Read radiative transfer rng state
  // ReadToString( dbfile, "rng_state_rad", proxy.rng_state_rad );
  // proxy.isReadRadiativeTransfer = true;

}

void AperData_IO::ReadParticleLoadMapFromOneFile( std::vector<size_t>& ptc_load_map, DBfile* dbfile, size_t (*load_calculator)(size_t N_e, size_t N_p, size_t N_i, size_t N_ph) ) const {
  DBtoc* dbtoc_snapshot = DBGetToc( dbfile );
  for ( int i = 0; i < dbtoc_snapshot->ndir; ++i ) {
    DBSetDir( dbfile, dbtoc_snapshot->dir_names[i] );
    // see if the current directory is one of an ensemble
    if ( DBInqVarExists(dbfile, "ens_label") ) {
      int ens_label = 0;
      DBReadVar( dbfile, "ens_label", &ens_label );

      // loop over all primary/ and replica/ directories for particle load
      DBtoc* dbtoc_ensemble = DBGetToc( dbfile );
      for ( int j = 0; j < dbtoc_ensemble->ndir; ++j ) {
        DBSetDir( dbfile, dbtoc_ensemble->dir_names[j] );

        auto f_readNum = [ dbfile ] ( std::string name ) {
          size_t num = 0;
          std::string varname = name + "/cell";
          if ( DBInqVarExists( dbfile, varname.c_str() ) ) {
            num = DBGetVarLength( dbfile, varname.c_str() );
          }
          return num;
        };
        size_t N_e = f_readNum( PtcType2Str(ParticleType::ELECTRON) );
        size_t N_p = f_readNum( PtcType2Str(ParticleType::POSITRON) );
        size_t N_i = f_readNum( PtcType2Str(ParticleType::ION) );
        size_t N_ph = f_readNum( PtcType2Str(ParticleType::PHOTON) );

        ptc_load_map[ens_label] += load_calculator( N_e, N_p, N_i, N_ph );

        DBSetDir( dbfile, ".." );
        dbtoc_ensemble = DBGetToc( dbfile );
      }
    }

    DBSetDir( dbfile, ".." );
    dbtoc_snapshot = DBGetToc( dbfile );
  }
}
