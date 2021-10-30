#pragma once
#include "dye/ensemble.hpp"
#include "field/sync.hpp"
#include "filesys/filesys.hpp"
#include "io/data_exporter.hpp"
#include "io/database_writer.hpp"

namespace io {
  namespace {
  template < typename Tcoarse, int DGrid, typename Tfine >
  apt::Grid<Tcoarse,DGrid> coarsen(const apt::Grid<Tfine,DGrid>& grid, int downsample_ratio) {
    apt::Grid<Tcoarse,DGrid> res;
    for ( int i = 0; i < DGrid; ++i ) {
      res[i] = apt::Grid1D<Tcoarse>( grid[i].lower(), grid[i].upper(), grid[i].dim() / downsample_ratio );
    }
    return res;
  }
  }


  template <typename RealDS, int DGrid, typename Real, template <typename> class S, typename RealJ>
  void DataExporter<RealDS, DGrid, Real, S, RealJ>::export_data(
      int timestep, Real dt, int num_files,
      const std::optional<mpi::CartComm> &cart_opt,
      const dye::Ensemble<DGrid> &ens,
      const apt::Grid<Real, DGrid> &grid, // local grid
      const field::Field<Real, 3, DGrid> &E,
      const field::Field<Real, 3, DGrid> &B,
      const field::Field<RealJ, 3, DGrid> &J, // J is Jmesh on a replica
      const particle::map<particle::array<Real, S>> &particles,
      const particle::map<particle::Properties> &properties) const {
    DatabaseWriter db_writer(cart_opt, m_data_dir, timestep, dt, num_files);
    ens.intra.barrier();

    db_writer.write_mesh<RealDS>( E.mesh(), m_is_collinear_mesh, m_mesh_ghost, m_downsample_ratio, grid );

    const auto grid_ds = coarsen<RealDS>( grid, m_downsample_ratio );

    field::Field<RealDS,3,DGrid> fds ( apt::make_range({}, apt::dims(grid_ds), m_mesh_ghost) );

    if ( cart_opt ) {
      for ( const auto& fe : m_fexps ) {
        fds = fe->action(m_downsample_ratio, cart_opt, grid, grid_ds, m_mesh_ghost, E, B, J );

        field::copy_sync_guard_cells( fds, *cart_opt );
        db_writer.write_var(fe->name(), "", fe->num_comps(), fds);
      }
    }

    for ( auto sp : particles ) {
      fds.reset();
      const auto& prop = properties[sp];

      for ( const auto& pe : m_pexps ) {
        fds = pe->action(m_downsample_ratio, grid, grid_ds, m_mesh_ghost, prop, particles[sp]);

        // FIXME reduce number of communication
        for (int i = 0; i < pe->num_comps(); ++i)
          ens.reduce_to_chief( mpi::by::SUM, fds[i].data().data(), fds[i].data().size() );
        if ( cart_opt ) {
          field::merge_sync_guard_cells( fds, *cart_opt );
          db_writer.write_var(pe->name(), prop.name, pe->num_comps(), fds);
        }
      }
    }
  }

}
