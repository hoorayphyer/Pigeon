#include "io/data_export.hpp"
#include "io/silo++.hpp"
#include "io/silo_optlist.hpp"
#include "parallel/mpi++.hpp"
#include "utility/filesys.hpp"
#include "field/communication.hpp"
#include "field/field_shape_interplay.hpp"

#include "parameters.hpp"
#include <unordered_map>

namespace io {
  std::string this_run_dir = "";

  void set_data_directory_for_this_run( std::string dir ) { this_run_dir = dir; }

  template < typename T, int DGrid >
  std::unordered_map<std::string, FieldBasedExportee<T,DGrid>*> fld_exportees;

  template < typename T, int DGrid >
  std::unordered_map<std::string, ParticleBasedExportee<T,DGrid>*> ptc_exportees;

  template < typename T, int DGrid >
  void register_exportee( std::string name, FieldBasedExportee<T,DGrid>* ptr ) {
    fld_exportees<T,DGrid>[name] = ptr;
  }

  template < typename T, int DGrid >
  void register_exportee( std::string name, ParticleBasedExportee<T,DGrid>* ptr ) {
    ptc_exportees<T,DGrid>[name] = ptr;
  }
}

namespace io {

  // TODO add one ghost zone in put_var
  template < typename T, int DGrid >
  void export_data( int timestep,
                    const Params& params,
                    const std::optional<mpi::Comm>& primary,
                    const std::optional<Ensemble<DGrid>>& ensemble,
                    const Mesh<T,DGrid>& export_mesh ) {

    if ( !ensemble ) return;

    constexpr int Padding = 6;
    char str_ts [Padding + 1];
    sprintf(str_ts, ("%0" + std::to_string(Padding)+ "d").c_str(), timestep);

    std::string dirname = filesys::append_slash(params.this_run_dir) + "data/timestep" + str_ts + "/";
    filesys::create_directories(dirname);

    silo::pmpio::file_t dbfile;
    expField<T,DGrid> exfd( data_grid );
    field::Field<T,1,DGrid> exfd( export_mesh );

    // TODO can we save just one mesh? In a separate file maybe?
    constexpr std::string meshname = "aperture_mesh";

    if ( primary ) {
      // TODO average to expf should factor in the scale functions, i.e. one should find the downsampled value by conserving the flux.

      dbfile = silo::pmpio::open<silo::Mode::Write>( dirname, *primary );

      {
        std::vector< std::vector<StorageType> > coords(DGrid);
        for ( int i = 0; i < Dim; ++i ) {
          auto& c = coords[i];
          auto dim = export_mesh.extent()[i];
          c.reserve(dim);
          c.resize(dim, static_cast<StorageType>(0));

          const auto& grid = export_mesh.bulk()[i];
          for (int j = 0; j < dim; ++j)
            c[j] = grid.absc( j, 0.5 );
        }

        silo::OptList optlist;
        optlist[DBOPT_TIME] = timestep * params.dt;
        optlist[DBOPT_CYCLE] = timestep;
        // optlist[DBOPT_BASEINDEX] = ; TODO is this needed?

        dbfile.put_mesh(meshname, coords, optlist);
      }

      for ( const auto [ name, exportee ] : fld_exportees<T,DGrid> ) {
        for ( const auto& I : apt::Block(export_mesh.bulk().extent()) ) {
          exfd[0]( I ) = (*exportee)( I, export_mesh );
        }
        sync_guard_cells_from_bulk( *primary, exfd );
        dbfile.put_var( name, meshname, exfd );
      }
    }

    {
      // TODO dbfile.put charge mass of the species in
      // TODO shapef should be with respect to datagrid. Rescaling?
      // FIXME TODO save gamma P density, which is total gamma / physical cell volume. Later this divided by number density gives us gamma P per unit particle. TODO: think this over.
      // TODO boundary conditions?? YES

      for ( const auto [ name, exportee ] : ptc_exportees<T,DGrid> ) {
        for ( auto& x : exfd[0].data() ) x = 0.0;
        auto beg = exportee->begin();
        auto end = exportee->end();
        for ( auto i = beg; i != end; exportee->next(i) )
          field::deposit( exfd, {exportee->val(i)}, exportee->loc(i), shapef );
        // TODO apply BC here

        ens.intra.reduce( ens.chief, exfd.data().data(), exfd.data.size() );
        if ( primary ) {
          merge_guard_cells_into_bulk( exfd, primary );
          sync_guard_cells_from_bulk( exfd, primary );
          dbfile.put_var( name, meshname, exfd );
        }
      }
    }

    silo::close( dbfile );

    {
      if ( !primary || primary->rank() != 0 ) return;

      // generate master file
      std::string masterfile = filesys::append_slash(params.this_run_dir) + "data/timestep" + str_ts + ".silo";

      auto dbmasterfile = silo::open<silo::Mode::Write>( masterfile );

      // dbmasterfile.put_multimesh ( timestep );

      // for (unsigned int i = 0; i < quadVars.size(); i++) {
      //   put_multivar( dbfileMaster, timestep, varname, vartype );
      // }

      // silo::close( dbfile );
    }

  }
}

