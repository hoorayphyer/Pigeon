#ifndef _IO_DATA_EXPORT_HPP_
#define _IO_DATA_EXPORT_HPP_

#include <string>
#include "io/exportee.hpp"

struct Params;
namespace mpi { struct Comm; }
namespace std {
  template <class T> class optional;
}

namespace io {
  void set_data_directory_for_this_run( std::string dir );

  template < typename T, int DGrid >
  void register_exportee( std::string name, FieldBasedExportee<T,DGrid>* );

  template < typename T, int DGrid >
  void register_exportee( std::string name, ParticleBasedExportee<T,DGrid>* );

  void export_data( int timestep,
                    const Params& params,
                    const std::optional<mpi::Comm>& primary,
                    const std::optional<mpi::Comm>& ensemble );
}

#endif
