#ifndef _AperData_IO_H_
#define _AperData_IO_H_
#include <string>
#include "AperData.h"

struct DBfile;
struct SaveSnapshotProxy;
class AperParams;
class MPIWorldCommunicator;
class MPICommunicator;
struct ProcessRoster;

struct AperData_IO {
  void WriteSnapshot( Scalar time, int timestep, const AperParams& params, const AperData& data, const MPICommunicator& comm, const ProcessRoster* roster, std::string snapshot_dir, const int num_pmpio_files ) const;

  void WriteTrackingData( int timestep, AperData& data, const MPIWorldCommunicator& world, const Grid& grid, std::string track_dir, const int num_pmpio_files ) const;

  void ReadSnapshot( AperData& data, const AperParams& params, const MPICommunicator& comm, std::string snapshot_dir ) const;

  int ReadTimestep( std::string snapshot_dir ) const;

  std::vector<int> ReadSizeMap( std::string snapshot_dir ) const;

  std::vector<size_t> ReadParticleLoadMap( std::string snapshot_dir, size_t (*load_calculator)(size_t N_e, size_t N_p, size_t N_i, size_t N_ph) ) const;

  friend class TestAperData_IO;

private:
  template <typename SinglePtc, typename Filter>
  void WriteParticle( DBfile* dbfile, const Particles<SinglePtc>& particles, Filter filter ) const;
  template <typename SinglePtc, typename Partitioner >
  void ReadParticle( DBfile* dbfile, Particles<SinglePtc>& particles, int ens_rank, const Partitioner& partitioner ) const;

  template <typename T, template <typename> class Field >
  void WriteField( DBfile* dbfile, const Field<T>& field, std::string fieldname ) const;
  template <typename T, template <typename> class Field >
  void ReadField( DBfile* dbfile, Field<T>& field, std::string fieldname ) const;

  void WriteTrackedParticles( DBfile* dbfile, const AperData& data, const Grid& grid, const std::string& track_dir, int world_rank ) const;
  void WriteTrackedPairCreation( DBfile* dbfile, PairCreationTracker& data, const std::string& track_dir, int world_rank ) const;

  void WriteExtra( DBfile* dbfile, const SaveSnapshotProxy& proxy ) const;
  void ReadExtra( DBfile* dbfile, SaveSnapshotProxy& proxy ) const;

  void ReadParticleLoadMapFromOneFile( std::vector<size_t>& ptc_load_map, DBfile* dbfile, size_t (*load_calculator)(size_t N_e, size_t N_p, size_t N_i, size_t N_ph) ) const;
};

#endif // ----- end of _AperData_IO_H_
