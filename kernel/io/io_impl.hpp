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
  template < int B, int E >
  constexpr int POW() {
    if constexpr ( E == 0 ) return 1;
    else return B * POW<B,E-1>();
  };

  constexpr char MeshExport[] = "PICMesh";

  template < int DGrid >
  constexpr auto ofs_export =
    []() {
      apt::array<field::offset_t, DGrid> res;
      apt::foreach<0,DGrid>( [](auto& a) { a = MIDWAY; }, res );
      return res;
    }();

  template < int DSRatio, int DGrid >
  constexpr apt::Index<DGrid> subext
  ( []() {
      apt::Index<DGrid> res;
      for ( int i = 0; i < DGrid; ++i ) res[i] = DSRatio;
      return res;}() );

  template < int DSRatio, int DGrid >
  constexpr auto Ifull( const apt::Index<DGrid>& Iexport, const apt::Index<DGrid>& Isub ) noexcept {
    apt::Index<DGrid> res;
    for ( int i = 0; i < DGrid; ++i ) res[i] = DSRatio * Iexport[i] + Isub[i];
    return res;
  }

  template < typename T, int DGrid >
  constexpr auto q_from_cell( const apt::Index<DGrid>& I, const apt::array<field::offset_t, DGrid>& ofs ) noexcept {
    apt::array<T,DGrid> res;
    for ( int i = 0; i < DGrid; ++i ) res[i] = I[i] + static_cast<T>(ofs[i]);
    return res;
  }

  // TODO This doesn't interpolate to zone center yet. Also may need intermediate sync_guards when interpolating to center. NOTE downsampling interpolation has its own shape function
  template < int DSRatio, typename RealExport, typename Real, int DGrid, int DField >
  void downsample ( field::Field<RealExport,1,DGrid>& io_field,
                    const field::Field<Real,DField,DGrid>& full_field, int comp ) {
    const auto& full_field_comp = full_field[comp];
    for ( const auto& I : apt::Block( io_field.mesh().bulk_dims() ) ) {
      Real f = 0.0;
      for ( const auto& Isub : apt::Block(subext<DSRatio, DGrid>) )
        f += full_field_comp( Ifull<DSRatio>(I,Isub) );
      io_field[0](I) = f / POW<DSRatio, DGrid>();
    }
  }

  template < typename Real, typename RealJ, int DGrid, typename Metric >
  void normalizeJ( field::Field<RealJ, 1, DGrid>& J, int comp, const mani::Grid<Real,DGrid>& grid ) {
    // NOTE J has DField 1 because it is meant to be an temporary field
    return;
  }

  template < int DSRatio,
             int DGrid,
             typename Real,
             template < typename > class PtcSpecs,
             typename ShapeF,
             typename RealJ,
             typename Metric >
  struct ExportEBJ {
    const mani::Grid<Real,DGrid>& grid; // local grid
    const field::Field<Real, 3, DGrid>& Efield;
    const field::Field<Real, 3, DGrid>& Bfield;
    const field::Field<RealJ, 3, DGrid>& Jfield;// J is Jmesh on a replica
    const particle::map<particle::array<Real,PtcSpecs>>& particles;

    template < typename RealExport, typename F_put_to_master >
    void operator() ( const silo::Pmpio& pmpio, field::Field<RealExport,1,DGrid>& io_field, const std::optional<mpi::CartComm>& cart_opt, const dye::Ensemble<DGrid>& ens, const F_put_to_master& put_to_master ) const {
      if ( !cart_opt ) return;

      field::Field<Real, 1, DGrid> tmp( Efield.mesh() );

      std::vector<int> dims(DGrid);
      for ( int i = 0; i < DGrid; ++i ) dims[i] = io_field.mesh().extent()[i];

      std::string varname;
      for ( int comp = 0; comp < 3; ++comp ) {
        {
          downsample<DSRatio>( io_field, Efield, comp );
          field::copy_sync_guard_cells( io_field, *cart_opt );
          varname = "E" + std::to_string(comp+1);
          pmpio([&](auto& dbfile){
                  dbfile.put_var( varname, MeshExport, io_field[0].data().data(), dims );
                });
          put_to_master( varname, cart_opt->size());

          downsample<DSRatio>( io_field, Bfield, comp );
          field::copy_sync_guard_cells( io_field, *cart_opt );
          varname = "B" + std::to_string(comp+1);
          pmpio([&](auto& dbfile){
                  dbfile.put_var( varname, MeshExport, io_field[0].data().data(), dims );
                });
          put_to_master( varname, cart_opt->size());
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
          downsample<DSRatio>( io_field, tmp, 0 );
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

  template < int DSRatio,
             int DGrid,
             typename Real,
             template < typename > class PtcSpecs,
             typename ShapeF,
             typename RealJ,
             typename Metric >
  struct ExportParticles {
  private:
    void fold_back_at_axis ( field::Field<Real,1,DGrid>& field ) const {
      constexpr int axis_dir = 1;
      constexpr auto PI = std::acos(-1.0l);
      bool is_at_axis_lower = std::abs( grid[axis_dir].lower() - 0.0 ) < grid[axis_dir].delta();
      bool is_at_axis_upper = std::abs( grid[axis_dir].upper() - PI ) < grid[axis_dir].delta();

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
    const mani::Grid<Real,DGrid>& grid; // local grid
    const field::Field<Real, 3, DGrid>& Efield;
    const field::Field<Real, 3, DGrid>& Bfield;
    const field::Field<RealJ, 3, DGrid>& Jfield;// J is Jmesh on a replica
    const particle::map<particle::array<Real,PtcSpecs>>& particles;

    template < typename RealExport, typename F_put_to_master >
    void operator() ( const silo::Pmpio& pmpio, field::Field<RealExport,1,DGrid>& io_field, const std::optional<mpi::CartComm>& cart_opt, const dye::Ensemble<DGrid>& ens, const F_put_to_master& put_to_master ) const {
      // TODO dbfile.put charge mass of the species in
      // TODO save gamma P density, which is total gamma / physical cell volume. Later this divided by number density gives us gamma P per unit particle. Think this over.
      field::Field<Real, 1, DGrid> tmp( Efield.mesh() );

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

          downsample<DSRatio>(io_field, tmp, 0);
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

          downsample<DSRatio>(io_field, tmp, 0);
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
             int DownsampleRatio,
             int DGrid,
             typename Real,
             template < typename > class PtcSpecs,
             typename RealJ
             >
  void export_data( std::string prefix, int timestep, Real dt, int num_files,
                    const std::optional<mpi::CartComm>& cart_opt,
                    const dye::Ensemble<DGrid>& ens,
                    const mani::Grid<Real,DGrid>& grid, // local grid
                    const field::Field<Real, 3, DGrid>& Efield,
                    const field::Field<Real, 3, DGrid>& Bfield,
                    const field::Field<RealJ, 3, DGrid>& Jfield,// J is Jmesh on a replica
                    const particle::map<particle::array<Real,PtcSpecs>>& particles
                    ) {
    constexpr int silo_mesh_ghost = 1;
    constexpr auto silo_mesh_type = silo::MeshType::Curv;

    char str_ts [10];
    sprintf(str_ts, "%06d\0", timestep);

    silo::file_t master; // only significant on world.rank() == 0
    if ( cart_opt && cart_opt->rank() == 0 ) {
      master = silo::open( prefix + "/timestep" + str_ts + ".silo", silo::Mode::Write );
    }

    prefix = prefix + "/data/timestep" + str_ts;

    silo::Pmpio pmpio; // only significant on carts
    if ( cart_opt ) {
      pmpio.filename = prefix + "/set" + std::to_string(cart_opt->rank() % num_files ) + ".silo";
      pmpio.dirname = "cart";
      for ( const auto& x : cart_opt->coords() ) {
        char tmp[10];
        sprintf(tmp, "_%03d", x );
        pmpio.dirname += tmp;
      }
      pmpio.comm = cart_opt->split( (cart_opt->rank()) % num_files );
      fs::mpido(*cart_opt, [&](){fs::create_directories(prefix);} );
    }

    ens.intra.barrier();

    // set up file_ns and block_ns. n below is thought of as the cartesian rank. Only significant on world.rank() == 0
    constexpr char delimiter = '|';
    // NOTE file namescheme uses relative path so that when the directory is moved, the data is still viewable
    const std::string file_ns = delimiter + std::string("data/timestep") + str_ts + "/set%d.silo" + delimiter + "n%" + std::to_string(num_files);

    auto block_ns_gen =
      [delimiter, &cart_opt]( std::string name ) -> std::string {
        static auto parts =
          [&]() -> apt::pair<std::string> {
            std::string part1 = delimiter + std::string("cart");
            auto [c, dims, p] = cart_opt->coords_dims_periodic();
            for ( int i = 0; i < dims.size(); ++i ) part1 += "_%03d";

            std::string part2 = "";
            // mpi uses row major numbering to map cartesian rank to coordinates, e.g. (0,0) -> 0, (0,1) -> 1, (1,0) -> 2, (1,1) -> 3
            std::vector<int> strides ( dims.size() + 1 ); // strides = ( DzDyDx, DzDy, Dz, 1 )
            strides.back() = 1;
            for ( int i = dims.size() - 1; i > -1; --i ) strides[i] = strides[i+1] * dims[i];

            for ( int i = 0; i < dims.size(); ++i ) {
              part2 += delimiter + std::string("(n%" + std::to_string(strides[i]) + ")/") + std::to_string(strides[i+1]);
            }
            return {part1, part2};
          } ();
        return parts[LFT] + "/" + name + parts[RGT];
      };

    auto bulk_dims_export = Efield.mesh().bulk_dims();
    for ( int i = 0; i < DGrid; ++i ) bulk_dims_export[i] /= DownsampleRatio;

    field::Field<RealExport,1,DGrid> field_export( { bulk_dims_export, silo_mesh_ghost } );

    // TODOL can we save just one mesh in a dedicated file and store a pointer to it in other files

    silo::OptList optlist;
    if ( cart_opt ) {
      optlist[silo::Opt::TIME] = timestep * dt;
      optlist[silo::Opt::CYCLE] = timestep;
    }

    auto put_to_master =
      [&cart_opt, &master, &file_ns, &block_ns_gen, &optlist]( std::string varname, int nblock ) {
        if ( cart_opt && cart_opt->rank() == 0 )
          master.put_multivar( varname, nblock, file_ns, block_ns_gen(varname), optlist );
      };

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
      auto optlist_mesh = optlist;
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
      //       c[j] = grid[i].absc( DownsampleRatio * j, ofs_export<DGrid>[i] );
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
          quadmesh_dims[i] = bulk_dims_export[i] + 1 + lo_ofs[i] + hi_ofs[i];

        RealExport* coords[DGrid];
        coords[0] = new RealExport [ quadmesh_dims[0] * quadmesh_dims[1] ];
        coords[1] = new RealExport [ quadmesh_dims[0] * quadmesh_dims[1] ];

        for ( int j = 0; j < quadmesh_dims[1]; ++j ) {
          for ( int i = 0; i < quadmesh_dims[0]; ++i ) {
            auto r = std::exp( grid[0].absc( DownsampleRatio * ( i - lo_ofs[0]) ) );
            auto theta = grid[1].absc( DownsampleRatio * (j - lo_ofs[1]) );
            coords[0][i + j * quadmesh_dims[0]] = r * std::sin(theta);
            coords[1][i + j * quadmesh_dims[0]] = r * std::cos(theta);
          }
        }

        pmpio( [&](auto& dbfile) {
                 DBPutQuadmesh(dbfile, MeshExport, NULL, coords, quadmesh_dims, DGrid, DB_FLOAT, DB_NONCOLLINEAR, optlist_mesh);
               } );

        delete [] coords[0];
        delete [] coords[1];

        if ( cart_opt->rank() == 0 )
          master.put_multimesh ( MeshExport, cart_opt->size(), file_ns, block_ns_gen(MeshExport), silo_mesh_type, optlist ); // NOTE use optlist here, not optlist_mesh
      }

      // TODOL
      // for ( const auto [ name, exportee ] : fld_exportees<T,DGrid> ) {
      //   for ( const auto& I : apt::Block(export_mesh.bulk().extent()) ) {
      //     exfd[0]( I ) = (*exportee)( I, export_mesh );
      //   }
      //   copy_sync_guard_cells( *primary, exfd );
      //   dbfile.put_var( name, meshname, exfd );
      // }

      ExportEBJ<DownsampleRatio, DGrid, Real, PtcSpecs, ShapeF, RealJ, Metric> exportEBJ{grid,Efield,Bfield,Jfield,particles};
      exportEBJ(pmpio, field_export, cart_opt, ens, put_to_master );
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

      ExportParticles<DownsampleRatio, DGrid, Real, PtcSpecs, ShapeF, RealJ, Metric> exportPtcs{grid,Efield,Bfield,Jfield,particles};
      exportPtcs(pmpio, field_export, cart_opt, ens, put_to_master );
    }

  }
}
