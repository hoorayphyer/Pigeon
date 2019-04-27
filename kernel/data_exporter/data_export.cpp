#include "data_export.hpp"
#include "io/silo++.hpp"
#include "parallel/mpi++.hpp"
#include "utility/filesys.hpp"
#include "field/communication.hpp"
#include "field/mesh_shape_interplay.hpp"

#include "dye/ensemble.hpp"

#include "gen.hpp"

#include <time.h>
#include <cstring>
namespace io {
  std::string this_run_dir;

  void init_this_run_dir( std::string prefix ) {
    util::fs::append_slash(prefix);

    char subDir[100] = {};
    for ( int i = 0; i < 100; ++i )
      subDir[i] = '\0';
    // use world root time to ensure uniqueness
    if ( mpi::world.rank() == 0 ) {
      char myTime[100] = {};
      time_t rawtime;
      struct tm* timeinfo;
      time (&rawtime);
      timeinfo = localtime(&rawtime);
      strftime(myTime, 100, "%Y%m%d-%H%M", timeinfo);
      snprintf(subDir, sizeof(subDir), "%s/", myTime);
    }
    mpi::world.broadcast( 0, subDir, 100 );

    this_run_dir = prefix + pic::project_name + "-" + subDir;

    if ( mpi::world.rank() == 0 ) {
      std::string local_data_dir = "Data/"; // local directory for storing data symlinks
      util::fs::create_directories(this_run_dir);
      util::fs::create_directories(local_data_dir);
      util::fs::create_directory_symlink(this_run_dir, local_data_dir + pic::project_name + "-" + subDir);
    }

    util::fs::append_slash(this_run_dir);
  }

  // TODO
  // void set_logger_dir( std::string logDir ) {
  //   // int rank = comm.world().rank();
  //   // Logger::thisRank = rank;
  //   // Logger::setActiveRank( rank );
  //   // if( Logger::isActiveRank && Logger::isLogToFile ) {
  //   //   FileSystem::create_directories(logDir);
  //   //   Logger::setLogFile( logDir + "rank_" + std::to_string( rank ) );
  //   // } else {
  //   //   // is not LogToFile, only let rank0 output to screen
  //   //   if ( rank != 0)
  //   //     Logger::setVerbosityLevel(-1);
  //   // }
  // }
}

namespace io {
  // template < typename T, int DGrid >
  // std::unordered_map<std::string, FieldBasedExportee<T,DGrid>*> fld_exportees;

  // template < typename T, int DGrid >
  // std::unordered_map<std::string, ParticleBasedExportee<T,DGrid>*> ptc_exportees;

  // template < typename T, int DGrid >
  // void register_exportee( std::string name, FieldBasedExportee<T,DGrid>* ptr ) {
  //   fld_exportees<T,DGrid>[name] = ptr;
  // }

  // template < typename T, int DGrid >
  // void register_exportee( std::string name, ParticleBasedExportee<T,DGrid>* ptr ) {
  //   ptc_exportees<T,DGrid>[name] = ptr;
  // }
}

namespace io {
  template < int B, int E >
  constexpr int POW() {
    if constexpr ( E == 0 ) return 1;
    else return B * POW<B,E-1>();
  };

  constexpr char MeshExport[] = "PICMesh";

  constexpr int ds_ratio = 2;
  using Tio = float;

  template < int DGrid >
  constexpr auto ofs_export =
    []() {
      apt::array<field::offset_t, DGrid> res;
      apt::foreach<0,DGrid>( [](auto& a) { a = MIDWAY; }, res );
      return res;
    }();

  template < int DGrid >
  constexpr apt::Index<DGrid> subext
  ( []() {
      apt::Index<DGrid> res;
      for ( int i = 0; i < DGrid; ++i ) res[i] = ds_ratio;
      return res;}() );

  template < int DGrid >
  constexpr auto Ifull( const apt::Index<DGrid>& Iexport, const apt::Index<DGrid>& Isub ) noexcept {
    apt::Index<DGrid> res;
    for ( int i = 0; i < DGrid; ++i ) res[i] = ds_ratio * Iexport[i] + Isub[i];
    return res;
  }

  template < typename T, int DGrid >
  constexpr auto q_from_cell( const apt::Index<DGrid>& I, const apt::array<field::offset_t, DGrid>& ofs ) noexcept {
    apt::array<T,DGrid> res;
    for ( int i = 0; i < DGrid; ++i ) res[i] = I[i] + static_cast<T>(ofs[i]);
    return res;
  }

  template < typename RealExport, typename Real, int DGrid, int DField, typename ShapeF >
  void downsample ( field::Field<RealExport,1,DGrid>& io_field,
                    const field::Field<Real,DField,DGrid>& full_field, int comp, const ShapeF& shapef ) {
    const auto& full_field_comp = full_field[comp];
    for ( const auto& I : apt::Block( io_field.mesh().bulk_dims()) ) {
      auto& f = io_field[0](I);
      f = 0.0;
      for ( const auto& Isub : apt::Block(subext<DGrid>) )
        f += field::interpolate( full_field_comp, q_from_cell<Real>( Ifull(I,Isub), ofs_export<DGrid> ), shapef );
      f /= POW<ds_ratio, DGrid>();
    }
  }

  template < typename Real, typename RealJ, int DGrid, typename Metric >
  void normalizeJ( field::Field<RealJ, 1, DGrid>& J, int comp, const knl::Grid<Real,DGrid>& grid ) {
    // NOTE J has DField 1 because it is meant to be an temporary field
    return;
  }

  template < int DGrid,
             typename Real,
             template < typename > class PtcSpecs,
             typename ShapeF,
             typename RealJ,
             typename Metric >
  struct ExportEBJ {
    const knl::Grid<Real,DGrid>& grid; // local grid
    const field::Field<Real, 3, DGrid>& Efield;
    const field::Field<Real, 3, DGrid>& Bfield;
    const field::Field<RealJ, 3, DGrid>& Jfield;// J is Jmesh on a replica
    const particle::map<particle::array<Real,PtcSpecs>>& particles;

    template < typename File_t, typename RealExport, typename F_put_to_master >
    void operator() ( File_t& dbfile, field::Field<RealExport,1,DGrid>& io_field, const std::optional<mpi::CartComm>& cart_opt, const dye::Ensemble<DGrid>& ens, const F_put_to_master& put_to_master ) const {
      if ( !cart_opt ) return;

      field::Field<Real, 1, DGrid> tmp( Efield.mesh() );

      std::vector<int> dims(DGrid);
      for ( int i = 0; i < DGrid; ++i ) dims[i] = io_field.mesh().extent()[i];


      std::string varname;
      for ( int comp = 0; comp < 3; ++comp ) {
        {
          downsample( io_field, Efield, comp, ShapeF() );
          field::sync_guard_cells_from_bulk( io_field, *cart_opt );
          varname = "E" + std::to_string(comp+1);
          dbfile.put_var( varname, MeshExport, io_field[0].data().data(), dims );
          if ( cart_opt->rank() == 0 ) put_to_master( varname, cart_opt->size());

          downsample( io_field, Bfield, comp, ShapeF() );
          field::sync_guard_cells_from_bulk( io_field, *cart_opt );
          varname = "B" + std::to_string(comp+1);
          dbfile.put_var( varname, MeshExport, io_field[0].data().data(), dims );
          if ( cart_opt->rank() == 0 ) put_to_master( varname, cart_opt->size());
        }

        {
          for ( const auto& I : apt::Block(tmp.mesh().bulk_dims()) ) {
            tmp[0](I) = Jfield[comp](I);
          }
          tmp.set_offset( 0, Jfield[comp].offset() );
          { // normalize J to orthonormal basis
            // define a function pointer.
            auto* h_func =
              [comp]() -> Real(*)(Real,Real,Real) {
                if ( comp < DGrid ) {
                  switch(comp) {
                  case 0: return Metric::template hh<0,Real>;
                  case 1: return Metric::template hh<1,Real>;
                  case 2: return Metric::template hh<2,Real>;
                  }
                } else {
                  return Metric::template hhh<Real>;
                }
              }();

            static_assert( DGrid == 2 );
            const auto& ofs = tmp[0].offset();
            apt::array<Real,DGrid> q {};
            for ( const auto& I : apt::Block(tmp.mesh().bulk_dims()) ) {
              // TODOL use generator
              for ( int i = 0; i < DGrid; ++i ) q[i] = grid[i].absc(I[i], ofs[i]);

              auto h = h_func(q[0], q[1], 0.0);
              if ( std::abs(h) > 1e-12 )
                tmp[0](I) /= h;
              else
                tmp[0](I) = 0.0;
            }
          }
          downsample( io_field, tmp, 0, ShapeF() );
          field::sync_guard_cells_from_bulk( io_field, *cart_opt );
          varname = "J"+std::to_string(comp+1);
          dbfile.put_var( varname, MeshExport, io_field[0].data().data(), dims );
          if ( cart_opt->rank() == 0 ) put_to_master( varname, cart_opt->size());
        }
      }

      // TODOL divE divB for tests. EdotB EdotJ for pulsar
    }
  };

  template < int DGrid,
             typename Real,
             template < typename > class PtcSpecs,
             typename ShapeF,
             typename RealJ,
             typename Metric >
  struct ExportParticles {
  private:
    template < typename RealExport >
    void fold_back_at_axis ( field::Field<RealExport,1,DGrid>& field ) const {
      constexpr int axis_dir = 1;
      bool is_at_axis_lower = std::abs( grid[axis_dir].lower() - 0.0 ) < grid[axis_dir].delta();
      bool is_at_axis_upper = std::abs( grid[axis_dir].upper() - pic::PI ) < grid[axis_dir].delta();

      const auto& mesh = field.mesh();
      const int guard = mesh.guard();

      // NOTE field is assumed to have all-MIDWAY offset
      if ( is_at_axis_lower ) {
        for ( const auto& trI : mesh.project(axis_dir, {}, mesh.extent() ) ) {
          for ( int n = 0; n < guard; ++n )
            field[0][ trI | n ] += field[0][ trI | -1 - n ];
        }
      }

      if ( is_at_axis_upper ) {
        const int dim = mesh.bulk_dim(axis_dir);
        for ( const auto& trI : mesh.project(axis_dir, {}, mesh.extent() ) ) {
          for ( int n = 0; n < guard; ++n )
            field[0][ trI | dim - 1 - n ] += field[0][ trI | dim + n ];
        }
      }
    }

  public:
    const knl::Grid<Real,DGrid>& grid; // local grid
    const field::Field<Real, 3, DGrid>& Efield;
    const field::Field<Real, 3, DGrid>& Bfield;
    const field::Field<RealJ, 3, DGrid>& Jfield;// J is Jmesh on a replica
    const particle::map<particle::array<Real,PtcSpecs>>& particles;

    template < typename File_t, typename RealExport, typename F_put_to_master >
    void operator() ( File_t& dbfile, field::Field<RealExport,1,DGrid>& io_field, const std::optional<mpi::CartComm>& cart_opt, const dye::Ensemble<DGrid>& ens, const F_put_to_master& put_to_master ) const {
      // TODO dbfile.put charge mass of the species in
      // FIXME TODO save gamma P density, which is total gamma / physical cell volume. Later this divided by number density gives us gamma P per unit particle. TODO: think this over.
      field::Field<Real, 1, DGrid> tmp( Efield.mesh() );

      std::vector<int> dims(DGrid);
      for ( int i = 0; i < DGrid; ++i ) dims[i] = io_field.mesh().extent()[i];

      // number density
      std::string varname;
      for ( const auto&[sp, ptcs] : particles ) {
        tmp.reset();
        const auto& prop = particle::properties[sp];
        for ( const auto& ptc : ptcs ) {
          field::deposit( tmp, {prop.charge_x}, field::to_standard( grid, ptc.q() ), ShapeF() );
        }
        fold_back_at_axis(tmp);
        downsample(io_field, tmp, 0, ShapeF());
        ens.reduce_to_chief( mpi::by::SUM, io_field[0].data().data(), io_field[0].data().size() );
        if ( cart_opt ) {
          field::sync_guard_cells_from_bulk( io_field, *cart_opt );
          varname = std::string("n") + particle::properties[sp].name;
          dbfile.put_var( varname, MeshExport, io_field[0].data().data(), dims );
          if ( cart_opt->rank() == 0 ) put_to_master( varname, cart_opt->size());
        }
      }
    }

  };

}


#include <silo.h> // for some DBOPTs
namespace io {

  template < typename RealExport,
             int DGrid,
             typename Real,
             template < typename > class PtcSpecs,
             typename ShapeF,
             typename RealJ,
             typename Metric >
  void export_data( int timestep, Real dt, int num_files,
                    const std::optional<mpi::CartComm>& cart_opt,
                    const dye::Ensemble<DGrid>& ens,
                    const knl::Grid<Real,DGrid>& grid, // local grid
                    const field::Field<Real, 3, DGrid>& Efield,
                    const field::Field<Real, 3, DGrid>& Bfield,
                    const field::Field<RealJ, 3, DGrid>& Jfield,// J is Jmesh on a replica
                    const particle::map<particle::array<Real,PtcSpecs>>& particles
                    ) {
    constexpr int export_guard = 1;  // add one ghost zone in put_var

    char str_ts [10];
    sprintf(str_ts, "%06d\0", timestep);

    const std::string prefix = this_run_dir + "data/timestep" + str_ts + "/";
    util::fs::create_directories(prefix);

    silo::pmpio::file_t dbfile;
    silo::file_t master;
    std::string file_ns;
    std::string block_ns;

    auto bulk_dims_export = Efield.mesh().bulk_dims();
    for ( int i = 0; i < DGrid; ++i ) bulk_dims_export[i] /= ds_ratio;
    field::Mesh<DGrid> mesh_export( bulk_dims_export, export_guard );

    field::Field<RealExport,1,DGrid> field_export( mesh_export );

    // TODOL can we save just one mesh in a dedicated file and store a pointer to it in other files

    // These are for push_mesh and put_multimesh
    silo::OptList optlist;
    std::vector<int> export_offset(DGrid);
    for ( auto& x : export_offset ) x = export_guard;

    optlist[DBOPT_TIME] = timestep * dt;
    optlist[DBOPT_CYCLE] = timestep;
    optlist[DBOPT_LO_OFFSET] = export_offset;
    optlist[DBOPT_HI_OFFSET] = export_offset;

    auto put_to_master =
      [&master, &file_ns, &block_ns, &optlist]( std::string varname, int nblock ) {
        // FIXME
        // master.put_multivar( varname, nblock, file_ns, block_ns, optlist );
      };

    if ( cart_opt ) {
      // TODOL average to expf should factor in the scale functions, i.e. one should find the downsampled value by conserving the flux.

      const int rank = cart_opt->rank();
      const auto coords = cart_opt->coords();
      std::string filename = prefix + "set" + std::to_string(rank % num_files)+".silo";
      std::string silo_dname = "cart";
      for ( const auto& x : coords ) {
        char tmp[10];
        sprintf(tmp, "_%03d", x );
        silo_dname += tmp;
      }

      dbfile = silo::pmpio::open<silo::Mode::Write>( filename, silo_dname, *cart_opt, num_files );

      if ( cart_opt->rank() == 0 ) {
        master = silo::open<silo::Mode::Write>( prefix + "../timestep" + str_ts + ".silo" );

        // set up file_ns and block_ns
        constexpr char delimiter = '|';
        file_ns = delimiter + prefix + "set%d.silo" + delimiter + "n%" + std::to_string(num_files);
        block_ns = delimiter + std::string("cart");
        auto [c, dims, p] = cart_opt -> coords_dims_periodic();
        for ( int i = 0; i < dims.size(); ++i ) block_ns += "_%03d";
        std::vector<int> strides ( dims.size() + 1 );
        strides[0] = 1;
        for ( int i = 0; i < dims.size(); ++i ) strides[i+1] = strides[i] * dims[i];
        for ( int i = 0; i < dims.size(); ++i ) {
          block_ns += delimiter + std::string("(n%" + std::to_string(strides[i+1]) + ")/") + std::to_string(strides[i]);
        }
      }
      {
        std::vector< std::vector<RealExport> > coords(DGrid);
        for ( int i = 0; i < DGrid; ++i ) {
          auto& c = coords[i];
          auto dim = mesh_export.extent()[i];
          c.reserve(dim);
          c.resize(dim, {});

          for ( int j = 0; j < dim; ++j )
            c[j] = grid[i].absc( ds_ratio * j, ofs_export<DGrid>[i] );
        }

        dbfile.put_mesh(MeshExport, coords, optlist);
        if ( cart_opt->rank() == 0 ) {
          master.put_multimesh ( MeshExport, cart_opt->size(), file_ns, block_ns, optlist );
        }
      }

      // TODOL
      // for ( const auto [ name, exportee ] : fld_exportees<T,DGrid> ) {
      //   for ( const auto& I : apt::Block(export_mesh.bulk().extent()) ) {
      //     exfd[0]( I ) = (*exportee)( I, export_mesh );
      //   }
      //   sync_guard_cells_from_bulk( *primary, exfd );
      //   dbfile.put_var( name, meshname, exfd );
      // }

      ExportEBJ<DGrid, Real, PtcSpecs, ShapeF, RealJ, Metric> exportEBJ{grid,Efield,Bfield,Jfield,particles};
      exportEBJ(dbfile, field_export, cart_opt, ens, put_to_master );
    }

    {
      // for ( const auto [ name, exportee ] : ptc_exportees<T,DGrid> ) {
      //   for ( auto& x : exfd[0].data() ) x = 0.0;
      //   auto beg = exportee->begin();
      //   auto end = exportee->end();
      //   for ( auto i = beg; i != end; exportee->next(i) )
      //     field::deposit( exfd, {exportee->val(i)}, exportee->loc(i), shapef );
      //   // apply BC here

      //   ens.intra.reduce( ens.chief, exfd.data().data(), exfd.data.size() );
      //   if ( primary ) {
      //     merge_guard_cells_into_bulk( exfd, primary );
      //     sync_guard_cells_from_bulk( exfd, primary );
      //     dbfile.put_var( name, meshname, exfd );
      //   }
      // }
      ExportParticles<DGrid, Real, PtcSpecs, ShapeF, RealJ, Metric> exportPtcs{grid,Efield,Bfield,Jfield,particles};
      exportPtcs(dbfile, field_export, cart_opt, ens, put_to_master );
    }

    silo::close( dbfile );
    silo::close( master );

  }
}

namespace io {
  using namespace pic;

  template
  void export_data<real_export_t, DGrid, real_t, particle::Specs, ShapeF, real_j_t, Metric>
  ( int timestep, real_t dt, int num_files,
    const std::optional<mpi::CartComm>& cart_opt,
    const dye::Ensemble<DGrid>& ens,
    const knl::Grid<real_t,DGrid>& grid, // local grid
    const field::Field<real_t, 3, DGrid>& Efield,
    const field::Field<real_t, 3, DGrid>& Bfield,
    const field::Field<real_j_t, 3, DGrid>& Jfield,// J is Jmesh on a replica
    const particle::map<particle::array<real_t,particle::Specs>>& particles
    );
}

