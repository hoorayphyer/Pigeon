#ifndef _IO_EXPORTER_HPP_
#define _IO_EXPORTER_HPP_

#include "io/exportee.hpp"
#include "io/saver.hpp"

namespace io {
  template < typename RDS, int DGrid, typename R, template < typename > class S, typename RJ >
  class DataExporter {
  private:
    int _ratio = 1;
    int _guard = 1;
    const std::optional<mpi::CartComm>& _cart_opt;
    const dye::Ensemble<DGrid>& _ens;

    using FexpT= FieldExportee<RDS, DGrid, R, RJ>;
    using PexpT = PtcExportee<RDS, DGrid, R, S>;
  public:

    DataExporter( int ratio, int guard, const std::optional<mpi::CartComm>& cart_opt, const dye::Ensemble<DGrid>& ens )
      : _ratio(ratio), _guard(guard), _cart_opt(cart_opt), _ens(ens) {}

    void execute( const DataSaver& saver,

                  const apt::Grid<R,DGrid>& grid,
                  const field::Field<R, 3, DGrid>& E,
                  const field::Field<R, 3, DGrid>& B,
                  const field::Field<RJ, 3, DGrid>& J,// J is Jmesh on a replica
                  const particle::map<particle::array<R,S>>& particles,
                  const particle::map<particle::Properties>& properties,
                  const std::vector<FexpT*>& fexps,
                  const std::vector<PexpT*>& Pexps
                  ) const;

  };
}

#endif
