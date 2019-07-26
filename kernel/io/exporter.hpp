#ifndef _IO_EXPORTER_HPP_
#define _IO_EXPORTER_HPP_

#include "io/exportee.hpp"
#include "io/saver.hpp"

namespace io {
  template < typename RealDS,
             int DGrid,
             typename Real,
             template < typename > class S,
             typename ShapeF,
             typename RealJ,
             typename Metric >
  class DataExporter {
  private:
    int _ratio = 1;
    int _guard = 1;
    const std::optional<mpi::CartComm>& _cart_opt;
    const dye::Ensemble<DGrid>& _ens;

  public:
    using FexpT= FieldExportee<RealDS, DGrid, Real, ShapeF, RealJ, Metric>;
    using PexpT = PtcExportee<RealDS, DGrid, Real, S, ShapeF>;

    DataExporter( int ratio, int guard, const std::optional<mpi::CartComm>& cart_opt, const dye::Ensemble<DGrid>& ens )
      : _ratio(ratio), _guard(guard), _cart_opt(cart_opt), _ens(ens) {}

    void execute( const DataSaver& saver,

                  const mani::Grid<Real,DGrid>& grid,
                  const field::Field<Real, 3, DGrid>& E,
                  const field::Field<Real, 3, DGrid>& B,
                  const field::Field<RealJ, 3, DGrid>& J,// J is Jmesh on a replica
                  const particle::map<particle::array<Real,S>>& particles,

                  const std::vector<FexpT*>& fexps,
                  const std::vector<PexpT*>& Pexps
                  ) const;

  };
}

#endif
