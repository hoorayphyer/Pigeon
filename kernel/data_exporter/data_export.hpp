#ifndef _IO_DATA_EXPORT_HPP_
#define _IO_DATA_EXPORT_HPP_

#include <string>
#include "kernel/grid.hpp"

#include "field/field.hpp"

#include "particle/map.hpp"
#include "particle/array.hpp"
// #include "io/exportee.hpp"

namespace mpi { struct CartComm; }
namespace std {
  template <class T> class optional;
}
namespace dye {
  template <int> struct Ensemble;
}

namespace io {
  // template < typename T, int DGrid >
  // void register_exportee( std::string name, FieldBasedExportee<T,DGrid>* );

  // template < typename T, int DGrid >
  // void register_exportee( std::string name, ParticleBasedExportee<T,DGrid>* );

  template < typename RealExport,
             int DGrid,
             typename Real,
             template < typename > class PtcSpecs,
             typename ShapeF,
             typename RealJ,
             typename Metric >
  void export_data( std::string this_run_dir, int timestep, Real dt, int num_files,
                    const std::optional<mpi::CartComm>& cart_opt,
                    const dye::Ensemble<DGrid>& ens,
                    const knl::Grid<Real,DGrid>& grid, // local grid
                    const field::Field<Real, 3, DGrid>& Efield,
                    const field::Field<Real, 3, DGrid>& Bfield,
                    const field::Field<RealJ, 3, DGrid>& Jfield,// J is Jmesh on a replica
                    const particle::map<particle::array<Real,PtcSpecs>>& particles
                    );
}

#endif
