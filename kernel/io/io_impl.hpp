#include "io/io.hpp"
#include "silopp/silo++.hpp"
#include "mpipp/mpi++.hpp"
#include "silopp/pmpio.hpp"
#include "filesys/filesys.hpp"
#include "field/sync.hpp"
#include "msh/mesh_shape_interplay.hpp"
#include "dye/ensemble.hpp"
#include "particle/properties.hpp"
#include <cmath>
#include "apt/numeric.hpp"
#include "manifold/curvilinear_impl.hpp"

#include "bc/axissymmetric.hpp"

#include "io/aux.hpp"

#include <silo.h> // for an ad hoc DBPutQuadmesh

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
  // local directory for storing data symlinks
#ifdef APPARENT_DATA_DIR
  std::string local_data_dir =
    []() {
      std::string str = APPARENT_DATA_DIR;
      fs::remove_slash(str);
      return str;
    }();
#else
  std::string local_data_dir = "Data";
#endif

  // TODOL what if prefix == local_data_dir??
  std::string init_this_run_dir( std::string prefix, std::string dirname ) {
    // use world root time to ensure uniqueness
    std::string this_run_dir;
    if ( mpi::world.rank() == 0 ) {
      prefix = fs::absolute(prefix);
      fs::remove_slash(prefix);
      fs::remove_slash(dirname);

      // in case of running too frequently within a minute, directories with postfixed numbers are created
      if ( fs::exists(prefix + "/" + dirname) ) {
        for ( int n = 1; ; ++n ) {
          if ( !fs::exists(prefix + "/" + dirname + "-" + std::to_string(n)) ) {
            dirname += "-" + std::to_string(n);
            break;
          }
        }
      }
      this_run_dir = prefix + "/" + dirname;

      fs::create_directories(this_run_dir);
      fs::create_directories(local_data_dir);
      fs::create_directory_symlink(this_run_dir, local_data_dir + "/" + dirname);
    }

    if ( mpi::world.size() > 1 ) {
      char buf[200];
      if ( mpi::world.rank() == 0 ) {
        for ( int i = 0; i < this_run_dir.size(); ++i )
          buf[i] = this_run_dir[i];
        buf[this_run_dir.size()] = '\0';
        mpi::world.broadcast(0, buf, 200);
      }
      else {
        mpi::world.broadcast(0, buf, 200);
        this_run_dir = {buf};
      }
    }

    return this_run_dir;

  }
}

namespace io {
  template < int DGrid,
             typename Real,
             template < typename > class S,
             typename ShapeF,
             typename RealJ,
             typename Metric >
  struct ExportEBJ {
    const mani::Grid<Real,DGrid>& grid; // local grid
    const field::Field<Real, 3, DGrid>& E;
    const field::Field<Real, 3, DGrid>& B;
    const field::Field<RealJ, 3, DGrid>& J;// J is Jmesh on a replica
    const particle::map<particle::array<Real,S>>& particles;

    template < typename RealExport, typename F_put_to_master >
    void operator() ( int DSRatio, const silo::Pmpio& pmpio, field::Field<RealExport,1,DGrid>& io_field, const std::optional<mpi::CartComm>& cart_opt, const dye::Ensemble<DGrid>& ens, const F_put_to_master& put_to_master ) const {
      if ( !cart_opt ) return;
      std::string MeshExport = "PICMesh";

      field::Field<Real, 1, DGrid> tmp( E.mesh() );

      std::vector<int> dims(DGrid);
      for ( int i = 0; i < DGrid; ++i ) dims[i] = io_field.mesh().extent()[i];

      std::string varname;
      for ( int comp = 0; comp < 3; ++comp ) {
        {
          downsample( DSRatio, io_field, E, comp );
          field::copy_sync_guard_cells( io_field, *cart_opt );
          varname = "E" + std::to_string(comp+1);
          pmpio([&](auto& dbfile){
                  dbfile.put_var( varname, MeshExport, io_field[0].data().data(), dims );
                });
          put_to_master( varname, cart_opt->size());

          downsample( DSRatio, io_field, B, comp );
          field::copy_sync_guard_cells( io_field, *cart_opt );
          varname = "B" + std::to_string(comp+1);
          pmpio([&](auto& dbfile){
                  dbfile.put_var( varname, MeshExport, io_field[0].data().data(), dims );
                });
          put_to_master( varname, cart_opt->size());
        }

        {
          for ( const auto& I : apt::Block(tmp.mesh().bulk_dims()) ) {
            tmp[0](I) = J[comp](I);
          }
          tmp.set_offset( 0, J[comp].offset() );
          { // normalize J to orthonormal basis
            // define a function pointer.
            Real(*h_func)(Real,Real,Real) = nullptr;
            switch(comp) {
            case 0: h_func = Metric::template hh<0,Real>; break;
            case 1: h_func = Metric::template hh<1,Real>; break;
            case 2: h_func = Metric::template hh<2,Real>; break;
            }

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
          downsample( DSRatio, io_field, tmp, 0 ); // 0 is to tmp
          field::copy_sync_guard_cells( io_field, *cart_opt );
          varname = "J"+std::to_string(comp+1);
          pmpio([&](auto& dbfile){
                  dbfile.put_var( varname, MeshExport, io_field[0].data().data(), dims );
                });
          put_to_master( varname, cart_opt->size());
        }
      }

      // TODOL divE divB for tests. EdotB EdotJ for pulsar
    }
  };

  template < int DGrid,
             typename Real,
             template < typename > class S,
             typename ShapeF,
             typename RealJ,
             typename Metric >
  struct ExportParticles {
  private:
    void fold_back_at_axis ( field::Field<Real,1,DGrid>& field ) const {
      // TODO these are duplicate
      constexpr int axis_dir = 1;
      constexpr auto PI = std::acos(-1.0l);
      bool is_lower = std::abs( grid[axis_dir].lower() - 0.0 ) < grid[axis_dir].delta();
      bool is_upper = std::abs( grid[axis_dir].upper() - PI ) < grid[axis_dir].delta();

      // NOTE field is assumed to have all-MIDWAY offset // TODO check this
      auto add_assign = []( Real& a, Real& b ) noexcept -> void {a += b; b = a;};
      bc::Axissymmetric<DGrid,Real,S,RealJ>::symmetrize( is_lower, is_upper, field[0], add_assign );
    }

  public:
    const mani::Grid<Real,DGrid>& grid; // local grid
    const field::Field<Real, 3, DGrid>& E;
    const field::Field<Real, 3, DGrid>& B;
    const field::Field<RealJ, 3, DGrid>& J;// J is Jmesh on a replica
    const particle::map<particle::array<Real,S>>& particles;

    template < typename RealExport, typename F_put_to_master >
    void operator() ( int DSRatio, const silo::Pmpio& pmpio, field::Field<RealExport,1,DGrid>& io_field, const std::optional<mpi::CartComm>& cart_opt, const dye::Ensemble<DGrid>& ens, const F_put_to_master& put_to_master ) const {
      // TODO dbfile.put charge mass of the species in
      // TODO save gamma P density, which is total gamma / physical cell volume. Later this divided by number density gives us gamma P per unit particle. Think this over.
      std::string MeshExport = "PICMesh";
      field::Field<Real, 1, DGrid> tmp( E.mesh() );

      std::vector<int> dims(DGrid);
      for ( int i = 0; i < DGrid; ++i ) dims[i] = io_field.mesh().extent()[i];

      std::string varname;
      // number density
      for ( const auto&[sp, ptcs] : particles ) {
        tmp.reset();
        const auto& prop = particle::properties.at(sp);
        for ( const auto& ptc : ptcs ) {
          if ( !ptc.is(particle::flag::exist) ) continue;
          msh::deposit( tmp, {1.0}, msh::to_standard( grid, ptc.q() ), ShapeF() );
        }
        // TODO to physical space
        ens.reduce_to_chief( mpi::by::SUM, tmp[0].data().data(), tmp[0].data().size() );
        if ( cart_opt ) {
          field::merge_sync_guard_cells( tmp, *cart_opt );
          fold_back_at_axis(tmp);

          downsample( DSRatio, io_field, tmp, 0);
          varname = std::string("n_") + particle::properties[sp].name;
          pmpio([&](auto& dbfile){
                  dbfile.put_var( varname, MeshExport, io_field[0].data().data(), dims );
                });
          put_to_master( varname, cart_opt->size());
        }
      }

      // Gamma density
      for ( const auto&[sp, ptcs] : particles ) {
        tmp.reset();
        const auto& prop = particle::properties.at(sp);
        for ( const auto& ptc : ptcs ) {
          if ( !ptc.is(particle::flag::exist) ) continue;
          msh::deposit( tmp, {std::sqrt( (prop.mass_x != 0) + apt::sqabs(ptc.p()) )}, msh::to_standard( grid, ptc.q() ), ShapeF() );
        }
        // TODO to physical space
        ens.reduce_to_chief( mpi::by::SUM, tmp[0].data().data(), tmp[0].data().size() );
        if ( cart_opt ) {
          field::merge_sync_guard_cells( tmp, *cart_opt );
          fold_back_at_axis(tmp);

          downsample(DSRatio, io_field, tmp, 0);
          varname = std::string("energy_") + particle::properties[sp].name;
          pmpio([&](auto& dbfile){
                  dbfile.put_var( varname, MeshExport, io_field[0].data().data(), dims );
                });
          put_to_master( varname, cart_opt->size());
        }
      }

    }
  };

}

namespace io {
  template < typename RealExport,
             typename Metric,
             typename ShapeF,
             int DGrid,
             typename Real,
             template < typename > class S,
             typename RealJ
             >
  void export_data( std::string prefix, int timestep, Real dt, int num_files, int downsample_ratio,
                    const std::optional<mpi::CartComm>& cart_opt,
                    const dye::Ensemble<DGrid>& ens,
                    const mani::Grid<Real,DGrid>& grid, // local grid
                    const field::Field<Real, 3, DGrid>& Efield,
                    const field::Field<Real, 3, DGrid>& Bfield,
                    const field::Field<RealJ, 3, DGrid>& Jfield,// J is Jmesh on a replica
                    const particle::map<particle::array<Real,S>>& particles
                    ) {
    constexpr int silo_mesh_ghost = 1;
    constexpr auto silo_mesh_type = silo::MeshType::Curv;

    char str_ts [10];
    sprintf(str_ts, "%06d\0", timestep);

    Downsampler<RealExport,DGrid> ds( downsample_ratio, Efield.mesh().bulk_dims(), silo_mesh_ghost, cart_opt );

    DataSaver saver(cart_opt);
    if ( cart_opt && cart_opt->rank() == 0 ) {
      saver.master = silo::open( prefix + "/timestep" + str_ts + ".silo", silo::Mode::Write );
      saver.set_namescheme( str_ts, num_files );
    }

    if ( cart_opt ) {
      saver.meshname = "PICMesh";

      saver.pmpio.filename = prefix + "/data/timestep" + str_ts + "/set" + std::to_string(cart_opt->rank() % num_files ) + ".silo";
      saver.pmpio.dirname = "cart";
      for ( const auto& x : cart_opt->coords() ) {
        char tmp[10];
        sprintf(tmp, "_%03d", x );
        saver.pmpio.dirname += tmp;
      }
      saver.pmpio.comm = cart_opt->split( (cart_opt->rank()) % num_files );
      fs::mpido(*cart_opt, [&](){fs::create_directories(prefix+ "/data/timestep" + str_ts);} );
    }

    ens.intra.barrier();

    // TODOL can we save just one mesh in a dedicated file and store a pointer to it in other files

    if ( cart_opt ) {
      saver.optlist[silo::Opt::TIME] = timestep * dt;
      saver.optlist[silo::Opt::CYCLE] = timestep;
    }

    if ( cart_opt ) {
      // set quadmesh ghost cells
      /* this is the most correct way to do ghost, but it needs mesh to support variable guards
      // std::vector<int> lo_ofs( DGrid, 0);
      // std::vector<int> hi_ofs( DGrid, 0);
      // {
      //   const auto& c = ens.cart_coords;
      //   for ( int i = 0; i < DGrid; ++i ) {
      //     if ( c[i] > 0 ) lo_ofs[i] = silo_mesh_ghost;
      //     if ( c[i] < ens.cart_dims[i] - 1 ) hi_ofs[i] = silo_mesh_ghost;
      //   }
      // }
      */
      std::vector<int> lo_ofs( DGrid, silo_mesh_ghost);
      std::vector<int> hi_ofs( DGrid, silo_mesh_ghost);
      auto optlist_mesh = saver.optlist;
      optlist_mesh[silo::Opt::LO_OFFSET] = lo_ofs;
      optlist_mesh[silo::Opt::HI_OFFSET] = hi_ofs;
      optlist_mesh[silo::Opt::BASEINDEX] = cart_opt -> coords(); // need this in rectilinear mesh

      // TODOL average to expf should factor in the scale functions, i.e. one should find the downsampled value by conserving the flux.

      // { // TODO uncomment after fixing the visit operator issue
      // TODO coord ghost cell issue
      //   std::vector< std::vector<RealExport> > coords(DGrid);
      //   for ( int i = 0; i < DGrid; ++i ) {
      //     auto& c = coords[i];
      //     auto dim = mesh_export.extent()[i];
      //     c.reserve(dim);
      //     c.resize(dim, {});

      //     for ( int j = 0; j < dim; ++j )
      //       c[j] = grid[i].absc( downsample_ratio * j, ofs_export<DGrid>[i] );
      //   }

      //   dbfile.put_mesh(MeshExport, coords, silo_mesh_type, optlist);
      //   if ( cart_opt->rank() == 0 ) {
      //     master.put_multimesh ( MeshExport, cart_opt->size(), file_ns, block_ns_gen(MeshExport), optlist_mesh );
      //   }
      // }

      { // TODO these are solely for LogSpherical2D. It is a hotfix on visit operators problems
        static_assert(DGrid == 2 );

        int quadmesh_dims[DGrid] = {};
        for ( int i = 0; i < DGrid; ++i )
          quadmesh_dims[i] = Efield.mesh().bulk_dim(i) / downsample_ratio + 1 + lo_ofs[i] + hi_ofs[i];

        RealExport* coords[DGrid];
        coords[0] = new RealExport [ quadmesh_dims[0] * quadmesh_dims[1] ];
        coords[1] = new RealExport [ quadmesh_dims[0] * quadmesh_dims[1] ];

        for ( int j = 0; j < quadmesh_dims[1]; ++j ) {
          for ( int i = 0; i < quadmesh_dims[0]; ++i ) {
            auto r = std::exp( grid[0].absc( downsample_ratio * ( i - lo_ofs[0]) ) );
            auto theta = grid[1].absc( downsample_ratio * (j - lo_ofs[1]) );
            coords[0][i + j * quadmesh_dims[0]] = r * std::sin(theta);
            coords[1][i + j * quadmesh_dims[0]] = r * std::cos(theta);
          }
        }

        saver.pmpio( [&](auto& dbfile) {
                       DBPutQuadmesh(dbfile, saver.meshname.c_str(), NULL, coords, quadmesh_dims, DGrid, DB_FLOAT, DB_NONCOLLINEAR, optlist_mesh);
               } );

        delete [] coords[0];
        delete [] coords[1];

        saver.PutMultimesh ( silo_mesh_type, saver.optlist ); // NOTE use saver.optlist here, not optlist_mesh
      }

      // ExportEBJ<DGrid, Real, S, ShapeF, RealJ, Metric> exportEBJ{grid,Efield,Bfield,Jfield,particles};
      // exportEBJ(downsample_ratio, pmpio, field_export, cart_opt, ens, put_to_master );
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
      //     merge_sync_guard_cells( exfd, primary );
      //     dbfile.put_var( name, meshname, exfd );
      //   }
      // }

      // ExportParticles<DGrid, Real, S, ShapeF, RealJ, Metric> exportPtcs{grid,Efield,Bfield,Jfield,particles};
      // exportPtcs(downsample_ratio, pmpio, field_export, cart_opt, ens, put_to_master );
    }

  }
}
