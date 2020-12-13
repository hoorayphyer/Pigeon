#include "io/io.hpp"
#include "silopp/silo++.hpp"
#include "mpipp/mpi++.hpp"
#include "silopp/pmpio.hpp"
#include "filesys/filesys.hpp"
#include "field/sync.hpp"
#include "dye/ensemble.hpp"
#include <cmath>

#include "io/exporter_impl.hpp"

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

// FIXME remove this
namespace io {
  bool is_collinear_mesh_silo = true;
  void set_is_collinear_mesh( bool x ) {
    is_collinear_mesh_silo = x;
  }
}

namespace io {
  template < typename RealDS,
             int DGrid,
             typename Real,
             template < typename > class S,
             typename RealJ
             >
  void export_data( std::string prefix, int timestep, Real dt, int num_files, int downsample_ratio,
                    const std::optional<mpi::CartComm>& cart_opt,
                    const dye::Ensemble<DGrid>& ens,
                    const apt::Grid<Real,DGrid>& grid, // local grid
                    const field::Field<Real, 3, DGrid>& E,
                    const field::Field<Real, 3, DGrid>& B,
                    const field::Field<RealJ, 3, DGrid>& J,// J is Jmesh on a replica
                    const particle::map<particle::array<Real,S>>& particles,
                    const particle::map<particle::Properties>& properties,
                    const std::vector<FieldExportee<RealDS, DGrid, Real, RealJ>*>& fexps,
                    const std::vector<PtcExportee<RealDS, DGrid, Real, S>*>& pexps
                    ) {
    using Exporter_t = DataExporter<RealDS, DGrid, Real, S, RealJ>;

    constexpr int silo_mesh_ghost = 1;
    auto silo_mesh_type = is_collinear_mesh_silo ? silo::MeshType::Rect : silo::MeshType::Curv;

    char str_ts [10];
    snprintf(str_ts, 10, "%07d", timestep);

    DataSaver saver(cart_opt);
    { // set up saver
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
        fs::mpido(*cart_opt, [&](){
                    fs::create_directories(prefix+ "/data/timestep" + str_ts);
                    saver.master.reset(new silo::file_t(
                        silo::open(prefix + "/timestep" + str_ts + ".silo",
                                   silo::Mode::Write)));
                    saver.set_namescheme(str_ts, num_files);
                  } );

        saver.optlist[silo::Opt::TIME] = timestep * dt;
        saver.optlist[silo::Opt::CYCLE] = timestep;
      }
    }

    ens.intra.barrier();

    Exporter_t exporter( downsample_ratio, silo_mesh_ghost, cart_opt, ens );

    // TODOL can we save just one mesh in a dedicated file and store a pointer to it in other files
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

      int quadmesh_dims[DGrid] = {};
      for ( int i = 0; i < DGrid; ++i )
        quadmesh_dims[i] = E.mesh().range(i).size() / downsample_ratio + 1 + lo_ofs[i] + hi_ofs[i]; // plus 1 to include the upper boundary

      if ( is_collinear_mesh_silo ) {
        std::vector< std::vector<RealDS> > coords(DGrid);
        for ( int i = 0; i < DGrid; ++i ) {
          auto& c = coords[i];
          auto dim = quadmesh_dims[i];
          c.reserve(dim);
          c.resize(dim, {});

          for ( int j = 0; j < dim; ++j )
            c[j] = grid[i].absc( downsample_ratio * j );
        }

        saver.pmpio( [&](auto& dbfile) {
                       dbfile.put_mesh(saver.meshname, coords, silo_mesh_type, optlist_mesh);
                     } );
      } else {
        // TODO these are solely for LogSpherical2D. It is a hotfix on visit operators problems
        static_assert(DGrid == 2 );

        RealDS* coords[DGrid];
        coords[0] = new RealDS [ quadmesh_dims[0] * quadmesh_dims[1] ];
        coords[1] = new RealDS [ quadmesh_dims[0] * quadmesh_dims[1] ];

        for ( int j = 0; j < quadmesh_dims[1]; ++j ) {
          for ( int i = 0; i < quadmesh_dims[0]; ++i ) {
            auto r = std::exp( grid[0].absc( downsample_ratio * ( i - lo_ofs[0]) ) );
            auto theta = grid[1].absc( downsample_ratio * (j - lo_ofs[1]) );
            coords[0][i + j * quadmesh_dims[0]] = r * std::sin(theta);
            coords[1][i + j * quadmesh_dims[0]] = r * std::cos(theta);
          }
        }

        saver.pmpio( [&](auto& dbfile) {
                       dbfile.put_mesh_noncollinear(saver.meshname, coords, quadmesh_dims, DGrid, optlist_mesh);
                     } );

        delete [] coords[0];
        delete [] coords[1];
      }
      saver.PutMultimesh ( silo_mesh_type, saver.optlist ); // NOTE use saver.optlist here, not optlist_mesh
    }

    exporter.execute( saver, grid, E, B, J, particles, properties, fexps, pexps );

  }
}
