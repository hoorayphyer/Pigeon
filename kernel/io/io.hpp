#ifndef _IO_DATA_EXPORT_HPP_
#define _IO_DATA_EXPORT_HPP_

#include <string>
#include "apt/grid.hpp"

#include "field/field.hpp"

#include "particle/map.hpp"
#include "particle/array.hpp"

#include "io/exportee.hpp"

namespace std {
  template <class T> class optional;
}
namespace dye {
  template <int> struct Ensemble;
}

namespace io {
  std::string init_this_run_dir( std::string prefix, std::string dirname );

  // TODO ad hoc
  void set_is_collinear_mesh( bool x );

  template < typename RDS,
             int DGrid,
             typename R,
             template < typename > class S,
             typename RJ
             >
  void export_data( std::string prefix, int timestep, R dt, int num_files, int downsample_ratio,
                    const std::optional<mpi::CartComm>& cart_opt,
                    const dye::Ensemble<DGrid>& ens,
                    const apt::Grid<R,DGrid>& grid, // local grid
                    const field::Field<R, 3, DGrid>& Efield,
                    const field::Field<R, 3, DGrid>& Bfield,
                    const field::Field<RJ, 3, DGrid>& Jfield,// J is Jmesh on a replica
                    const particle::map<particle::array<R,S>>& particles,
                    const particle::map<particle::Properties>& properties,
                    const std::vector<FieldExportee<RDS, DGrid, R, RJ>*>& fexps,
                    const std::vector<PtcExportee<RDS, DGrid, R, S>*>& pexps
                    );
}

#endif
