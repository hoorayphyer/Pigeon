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

#include "io/exporter_impl.hpp"

#include <silo.h> // for an ad hoc DBPutQuadmesh

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
  template < int F, typename Real, int DGrid, typename RealJ >
  apt::array<Real,3> field_self ( apt::Index<DGrid> I,
                                  const mani::Grid<Real,DGrid>& grid,
                                  const field::Field<Real, 3, DGrid>& E,
                                  const field::Field<Real, 3, DGrid>& B,
                                  const field::Field<RealJ, 3, DGrid>& J) {
    if constexpr ( F == 0 ) {
        return { E[0](I), E[1](I), E[2](I) };
      } else if ( F == 1 ) {
      return { B[0](I), B[1](I), B[2](I) };
    } else if ( F == 2 ) {
      return { J[0](I), J[1](I), J[2](I) };
    } else {
      static_assert(F < 3);
    }
  }

  // TODO add interpolate to center for E, B, J

  template < typename RealDS, int DGrid, typename Metric >
  void divide_flux_by_area ( field::Field<RealDS,3,DGrid>& fds, const mani::Grid<RealDS,DGrid>& grid, int num_comps = 3  ) {
    // define a function pointer.
    RealDS(*h_func)(RealDS,RealDS,RealDS) = nullptr;
    apt::array<RealDS,3> q {}; // TODO 3 is hard coded

    for ( int comp = 0; comp < num_comps; ++comp ) {
      const auto& ofs = fds[comp].offset();
      switch(comp) {
      case 0: h_func = Metric::template hh<0,RealDS>; break;
      case 1: h_func = Metric::template hh<1,RealDS>; break;
      case 2: h_func = Metric::template hh<2,RealDS>; break;
      }

      for ( const auto& I : apt::Block(fds.mesh().bulk_dims()) ) {

        for ( int i = 0; i < DGrid; ++i ) q[i] = grid[i].absc(I[i], ofs[i]);

        auto h = h_func(q[0], q[1], q[2]);
        if ( std::abs(h) > 1e-12 )
          fds[comp](I) /= h;
        else
          fds[comp](I) = 0.0;
      }
    }
  }

  template < int DGrid, typename Real, template < typename > class S, typename RealJ >
  void fold_back_at_axis ( field::Field<Real,3,DGrid>& field, const mani::Grid<Real,DGrid>& grid, int num_comps ) {
    // NOTE field is assumed to have all-MIDWAY offset
    for ( int i = 0; i < num_comps; ++i ) {
      for ( int dim = 0; dim < DGrid; ++dim )
        assert( field[i].offset()[dim] == MIDWAY );
    }

    // TODO these are duplicate
    constexpr int axis_dir = 1;
    constexpr auto PI = std::acos(-1.0l);
    bool is_lower = std::abs( grid[axis_dir].lower() - 0.0 ) < grid[axis_dir].delta();
    bool is_upper = std::abs( grid[axis_dir].upper() - PI ) < grid[axis_dir].delta();

    auto add_assign = []( Real& a, Real& b ) noexcept -> void {a += b; b = a;}; // TODO only add_assign
    for ( int i = 0; i < num_comps; ++i )
      bc::Axissymmetric<DGrid,Real,S,RealJ>::symmetrize( is_lower, is_upper, field[i], add_assign );
  }
}

namespace io {
  template < typename Real, template < typename > class S >
  apt::array<Real,3> ptc_num ( const particle::Properties& prop, const typename particle::array<Real,S>::const_particle_type& ptc ) {
    return { 1.0, 0.0, 0.0 };
  }

  template < typename Real, template < typename > class S >
  apt::array<Real,3> ptc_energy ( const particle::Properties& prop, const typename particle::array<Real,S>::const_particle_type& ptc ) {
    return { std::sqrt( (prop.mass_x != 0) + apt::sqabs(ptc.p()) ), 0.0, 0.0 };
  }

  template < typename Real, template < typename > class S >
  apt::array<Real,3> ptc_momentum ( const particle::Properties& prop, const typename particle::array<Real,S>::const_particle_type& ptc ) {
    return { ptc.p()[0], ptc.p()[1], ptc.p()[2] };
  }
}

namespace io {
  template < typename RealDS,
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
                    const field::Field<Real, 3, DGrid>& E,
                    const field::Field<Real, 3, DGrid>& B,
                    const field::Field<RealJ, 3, DGrid>& J,// J is Jmesh on a replica
                    const particle::map<particle::array<Real,S>>& particles
                    ) {
    using Exporter_t = DataExporter<RealDS, DGrid, Real, S, ShapeF, RealJ, Metric>;

    constexpr int silo_mesh_ghost = 1;
    constexpr auto silo_mesh_type = silo::MeshType::Curv;

    char str_ts [10];
    sprintf(str_ts, "%06d\0", timestep);

    DataSaver saver(cart_opt);
    { // set up saver
      if ( cart_opt && cart_opt->rank() == 0 ) {
        saver.master.reset( new silo::file_t( silo::open( prefix + "/timestep" + str_ts + ".silo", silo::Mode::Write ) ) );
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

      // { // FIXME uncomment after fixing the visit operator issue
      //   std::vector< std::vector<RealDS> > coords(DGrid);
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
          quadmesh_dims[i] = E.mesh().bulk_dim(i) / downsample_ratio + 1 + lo_ofs[i] + hi_ofs[i];

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
                       DBPutQuadmesh(dbfile, saver.meshname.c_str(), NULL, coords, quadmesh_dims, DGrid, DB_FLOAT, DB_NONCOLLINEAR, optlist_mesh);
               } );

        delete [] coords[0];
        delete [] coords[1];

        saver.PutMultimesh ( silo_mesh_type, saver.optlist ); // NOTE use saver.optlist here, not optlist_mesh
      }
    }

    std::vector<typename Exporter_t::FexpT*> fexps;
    std::vector<typename Exporter_t::PexpT*> pexps;

    {
      using FA = FieldAction<RealDS,DGrid,Real,RealJ,Metric>;

      fexps.push_back( new FA ( "E", 3,
                                field_self<0,Real, DGrid, RealJ>,
                                nullptr
                                ) );
      fexps.push_back( new FA ( "B", 3,
                                field_self<1,Real, DGrid, RealJ>,
                                nullptr
                                ) );
      fexps.push_back( new FA ( "J", 3,
                                field_self<2,Real, DGrid, RealJ>,
                                divide_flux_by_area<RealDS, DGrid, Metric>
                                ) );
    }

    {
      using PA = PtcAction<RealDS,DGrid,Real,S,ShapeF>;
      pexps.push_back( new PA ("Num", 1,
                               ptc_num<Real,S>,
                               fold_back_at_axis< DGrid, Real, S, RealJ >
                               ) );

      pexps.push_back( new PA ("E", 1,
                               ptc_energy<Real,S>,
                               fold_back_at_axis< DGrid, Real, S, RealJ >
                               ) );

      pexps.push_back( new PA ("P", 3,
                               ptc_energy<Real,S>,
                               fold_back_at_axis< DGrid, Real, S, RealJ >
                               ) );
    }


    exporter.execute( saver, grid, E, B, J, particles, fexps, pexps );

    for ( auto ptr : fexps ) delete ptr;
    for ( auto ptr : pexps ) delete ptr;

  }
}
