#include "ckpt/checkpoint.hpp"
#include "mpipp/mpi++.hpp"
#include "silopp/silo++.hpp"
#include "silopp/pmpio.hpp"
#include "filesys/filesys.hpp"
#include "particle/properties.hpp"
#include "field/field.hpp"
#include <cassert>
#include "dye/dynamic_balance.hpp"
#include "dye/scatter_load.hpp"

#ifdef PIC_DEBUG
#include "debug/debugger.hpp"
#include "logger/logger.hpp"
#endif

// auxilliary
namespace ckpt {
  template < typename T, int DGrid >
  struct FieldCkpt {
    template < int DField >
    void save( silo::file_t& dbfile, const std::string& name, const field::Field<T,DField,DGrid>& f ) {
      for ( int comp = 0; comp < DField; ++comp ) {
        auto n = name + ( DField == 1 ? "" : std::to_string(comp+1) );
        dbfile.write( n, f._comps[comp] );
      }
    }

    template < int DField >
    void load( silo::file_t& dbfile, const std::string& name, field::Field<T,DField,DGrid>& f ) {
      for ( int comp = 0; comp < DField; ++comp ) {
        auto n = name + ( DField == 1 ? "" : std::to_string(comp+1) );
#ifdef PIC_DEBUG
        lgr::file << "LDCKPT Field " << n << ", length of ckpt data = " << dbfile.var_length(n) << ", length of container = " << f._comps[comp].size() << std::endl;
#endif
        dbfile.read( n, f._comps[comp].data() );
#ifdef PIC_DEBUG
        lgr::file << silo::errmsg() << std::endl;
#endif
      }
    }
  };

  template < typename T, template < typename > class S >
  struct ParticleArrayCkpt {
    void save( silo::file_t& dbfile, const particle::Properties& prop, const particle::array<T,S>& ptcs ) {
      dbfile.mkcd(prop.name);
      particle::load_t num = ptcs.size();
      // NOTE: explicitly write load 1) to signal that this species in principle can exist and 2) because in Silo array length has type int so inferring from it might in some extreme case is wrong and hard to debug
      dbfile.write("load", num);
      if ( num != 0 ) {
        constexpr auto DPtc = S<T>::Dim;
        for ( int i = 0; i < DPtc; ++i ) {
          dbfile.write("q"+std::to_string(i+1), ptcs._q[i].data(), {ptcs._q[i].size()} );
          dbfile.write("p"+std::to_string(i+1), ptcs._p[i].data(), {ptcs._p[i].size()} );
        }
        dbfile.write( "frac", ptcs._frac.data(), {ptcs._frac.size()} );
        dbfile.write( "state", ptcs._state.data(), {ptcs._state.size()} );
      }
      dbfile.cd("..");
    }

    void load( silo::file_t& dbfile, const particle::Properties& prop, particle::array<T,S>& ptcs, int receiver_idx, int num_receivers ) {
      dbfile.cd( prop.name );
      auto N = dbfile.read1<particle::load_t>("load");
#ifdef PIC_DEBUG
      lgr::file << "LDCKPT SPECIES " << prop.name << ", Load = " << N << std::endl;
#endif

      if ( 0 != N ) {
        auto[ from, num ] = dye::scatter_load(N, receiver_idx, num_receivers);
        auto size_old = ptcs.size();
        ptcs.resize(size_old + num);

        constexpr auto DPtc = S<T>::Dim;
        apt::array<int, 3> slice { from, num, 1 };
#ifdef PIC_DEBUG
        lgr::file << "LDCKPT Slice from = " << from << ", length = " << num << ", stride = " << 1 << std::endl;
#endif
        {
          for ( int i = 0; i < DPtc; ++i ) {
            dbfile.readslice("q"+std::to_string(i+1), {slice}, ptcs._q[i].data() + size_old );
#ifdef PIC_DEBUG
            lgr::file << silo::errmsg() << std::endl;
#endif
            dbfile.readslice("p"+std::to_string(i+1), {slice}, ptcs._p[i].data() + size_old );
#ifdef PIC_DEBUG
            lgr::file << silo::errmsg() << std::endl;
#endif
          }
          dbfile.readslice("frac", {slice}, ptcs._frac.data() + size_old );
          dbfile.readslice("state", {slice}, ptcs._state.data() + size_old );
#ifdef PIC_DEBUG
          lgr::file << silo::errmsg() << std::endl;
#endif
        }
      }
      dbfile.cd("..");
    }

    inline const auto& q( const particle::array<T,S>& ptcs, int i ) const noexcept { return ptcs._q[i]; }
    inline const auto& p( const particle::array<T,S>& ptcs, int i ) const noexcept { return ptcs._p[i]; }
    inline const auto& frac( const particle::array<T,S>& ptcs ) const noexcept { return ptcs._frac; }
    inline const auto& state( const particle::array<T,S>& ptcs ) const noexcept { return ptcs._state; }

  };
}

namespace ckpt {
  // TODO save random seed
  template < int DGrid, typename R, template < typename > class S >
  std::string save_checkpoint( std::string prefix, const int num_parts,
                               const std::optional<dye::Ensemble<DGrid>>& ens_opt,
                               int timestep,
                               const field::Field<R, 3, DGrid>& E,
                               const field::Field<R, 3, DGrid>& B,
                               const particle::map<particle::array<R,S>>& particles,
                               const particle::map<particle::Properties>& properties,
                               const particle::map<R>& N_scat
                               ) {
    bool is_idle = !ens_opt;
    int key = 0;
    if ( !is_idle ) key = ens_opt -> label();
    // use ensemble label to put processes from same ensemble together
    auto active = mpi::world.split( {is_idle}, key );
    if ( is_idle ) return "";
    const auto& ens = *ens_opt;

    char str_ts [10];
    sprintf(str_ts, "%06d\0", timestep);

    prefix = prefix + "/checkpoints/timestep" + str_ts;

    silo::Pmpio pmpio;

    int num_per_part = active->size() / num_parts + ( active->size() % num_parts != 0 );

    {
      pmpio.filename = prefix + "/part" + std::to_string( active->rank() / num_per_part ) + ".silo";
      pmpio.dirname = "/ensemble" + std::to_string(ens.label());

      pmpio.comm = active -> split( active->rank() / num_per_part );
      fs::mpido(*active, [&](){fs::create_directories(prefix);} );
    }

    active->barrier();

    std::vector<particle::load_t> loads;
    loads.reserve(particles.size());
    for ( auto sp : particles )
      loads.push_back(particles[sp].size());

    ens.reduce_to_chief( mpi::by::SUM, loads.data(), loads.size() );

    pmpio([&]( auto& dbfile )
          {
            // write global data
            if ( !dbfile.exists("/timestep") ) {
              dbfile.write( "/cartesian_topology", reinterpret_cast<const int*>(ens.cart_topos.begin()), {DGrid} );

              using T = std::underlying_type_t<particle::species>;
              std::vector<T> sps;
              for ( auto sp : particles )
                sps.push_back( static_cast<T>(sp) );
              dbfile.write( "/species", sps  );

              dbfile.write( "/timestep", timestep );
            }

            { // write ensemble-wise data
              if ( !dbfile.exists("label") ) {
                dbfile.write( "cartesian_coordinates", ens.cart_coords.begin(), {DGrid} );
                dbfile.write( "label", ens.label() );
              }

              FieldCkpt<R,DGrid> ckpt;
              if ( ens.intra.rank() == ens.chief ) {
                int idx = 0;
                for ( auto sp : particles )
                  dbfile.write( properties[sp].name + "_load", loads[idx++] );

                ckpt.save(dbfile, "E", E);
                ckpt.save(dbfile, "B", B);

                { // write N_scat
                  std::vector<int> buffer;
                  for ( auto sp : N_scat ) buffer.push_back(static_cast<int>(sp));
                  dbfile.write( "N_scat_sp", buffer );
                  dbfile.write( "N_scat_data", N_scat.data() );
                }
              }
            }

            dbfile.mkcd( "rank" + std::to_string(ens.intra.rank()) );
            { // write process specific data
              int r = ens.intra.rank();
              dbfile.write( "r", r );
              ParticleArrayCkpt<R, S> ckpt;
              for ( auto sp : particles ) {
                ckpt.save( dbfile, properties[sp], particles[sp] );
              }
            }
            dbfile.cd("..");
          }
      );

    return prefix;

  }

  template < int DGrid, typename R, template < typename > class S >
  int load_checkpoint( std::string dir,
                       std::optional<dye::Ensemble<DGrid>>& ens_opt,
                       const std::optional<mpi::CartComm>& cart_opt,
                       field::Field<R, 3, DGrid>& E,
                       field::Field<R, 3, DGrid>& B,
                       particle::map<particle::array<R,S>>& particles,
                       const particle::map<particle::Properties>& properties,
                       particle::map<R>& N_scat,
                       int target_load
                       ) {
    int checkpoint_ts = 0;
    using T = std::underlying_type_t<particle::species>;
    int num_sps = 0;
    std::vector<T> sps;
    int num_ens = 0;
    particle::load_t myload = 0;

    // rank0 read data
    if ( mpi::world.rank() == 0 ) {
      auto dir_itr = fs::directory_iterator(dir);
      auto dbfile = silo::open(*dir_itr, silo::Mode::Read);
      assert( dbfile.var_exists("/timestep") );
      assert( dbfile.var_exists("/cartesian_topology") );
      assert( dbfile.var_exists("/species") );

      dbfile.read("/timestep", &checkpoint_ts);
#ifdef PIC_DEBUG
      lgr::file << "LDCKPT timestep = " << checkpoint_ts << std::endl;
      lgr::file << silo::errmsg() << std::endl;
#endif
      {
        int ndims = dbfile.var_length("/cartesian_topology");
        std::vector<int> topos(ndims);
        dbfile.read("/cartesian_topology", topos.data());
#ifdef PIC_DEBUG
        lgr::file << "LDCKPT cartesian_topology = (";
        for ( auto x : topos ) lgr::file << x << ",";
        lgr::file << ")" << std::endl;
        lgr::file << silo::errmsg() << std::endl;
#endif
        num_ens = 1;
        for ( auto x : topos ) num_ens *= std::abs(x);
      }
      {
        num_sps = dbfile.var_length("/species");
        sps.resize( num_sps );
        dbfile.read("/species", sps.data() );
#ifdef PIC_DEBUG
        lgr::file << "LDCKPT species = (";
        for ( auto x : sps ) lgr::file << x << ",";
        lgr::file << ")" << std::endl;
        lgr::file << silo::errmsg() << std::endl;
#endif
      }
    }
    // rank0 broadcast data
    {
      int buf[3] = {checkpoint_ts, num_sps, num_ens };
      mpi::world.broadcast( 0, buf, 3 );
      if ( mpi::world.rank() != 0 ) {
        checkpoint_ts = buf[0];
        num_sps = buf[1];
        num_ens = buf[2];
        sps.resize(num_sps);
      }
      mpi::world.broadcast( 0, sps.data(), sps.size() );
#ifdef PIC_DEBUG
      if ( mpi::world.rank() != 0 ) {
        lgr::file << "LDCKPT timestep = " << checkpoint_ts << std::endl;

        lgr::file << "LDCKPT num_ens = " << num_ens << std::endl;

        lgr::file << "LDCKPT species = (";
        for ( auto x : sps ) lgr::file << x << ",";
        lgr::file << ")" << std::endl;
      }
#endif
    }
    // each primary figure out its load
    if ( cart_opt ) {
      assert( ens_opt );
      const int mylabel = ens_opt->label();
      for ( auto f : fs::directory_iterator(dir) ) {
        bool is_found = false;
        auto sf = silo::open(f, silo::Mode::Read);
        for ( const auto& dname : sf.toc_dir() ) {
          if ( sf.var_exists( dname + "/rank0" ) ) {
            int l = sf.read1<int>(dname+"/label");
            if ( mylabel == l ) {
              sf.cd(dname);
              for ( int i = 0; i < num_sps; ++i ) {
                auto sp = static_cast<particle::species>(sps[i]);
                std::string entry = properties[sp].name + "_load";
                assert( sf.var_exists(entry) );
                myload += sf.read1<particle::load_t>(entry);
              }
#ifdef PIC_DEBUG
              lgr::file << "LDCKPT Label = " << mylabel << ", Load = " << myload << std::endl;
              lgr::file << silo::errmsg() << std::endl;
#endif
              is_found = true;
              sf.cd("..");
              break;
            }
          }
        }
        if ( is_found ) break;
      }
    }

    // dynamic assign computational resources
    {
      auto ens_opt_new = dye::deploy(myload, ens_opt, cart_opt, target_load);
      ens_opt.swap(ens_opt_new);
    }

    // by this step, all processes should have ens_opt set properly. The real loading begins now
    if ( !ens_opt ) return checkpoint_ts;

    const int mylabel = ens_opt->label();
    const int myrank = ens_opt->intra.rank();
    const auto ens_size = ens_opt->intra.size();

    FieldCkpt<R,DGrid> f_ckpt;
    ParticleArrayCkpt<R, S> p_ckpt;
    for ( auto f : fs::directory_iterator(dir) ) {
#ifdef PIC_DEBUG
      lgr::file << "Reading file " << f << std::endl;
#endif
      auto sf = silo::open(f, silo::Mode::Read);
      for ( const auto& dname : sf.toc_dir() ) {
        if ( dname.find("ensemble") != 0 ) continue;
        int l = sf.read1<int>(dname+"/label");
        if ( l != mylabel ) continue;
#ifdef PIC_DEBUG
        lgr::file << "Reading EB from " << dname << std::endl;
#endif
        sf.cd(dname);
        if ( 0 == myrank && sf.var_exists("rank0") ) { // NOTE: ranks of one ensemble may exist across files
          f_ckpt.load( sf, "E", E );
          f_ckpt.load( sf, "B", B );
          {
            assert( sf.var_exists("N_scat_sp") );
            assert( sf.var_exists("N_scat_data") );

            const auto buf_sp = sf.read1d<int>("N_scat_sp");
#ifdef PIC_DEBUG
            lgr::file << "LDCKPT N_scat_sp, size = " << sf.var_length("N_scat_sp") << ", data = (";
            for ( auto x : buf_sp ) lgr::file << x << ", ";
            lgr::file << ")" << std::endl << silo::errmsg() << std::endl;
#endif

            const auto buf_data = sf.read1d<double>("N_scat_data");
#ifdef PIC_DEBUG
            lgr::file << "LDCKPT N_scat_data, size = " << sf.var_length("N_scat_data") << ", data = (";
            for ( auto x : buf_data ) lgr::file << x << ", ";
            lgr::file << ")" << std::endl << silo::errmsg() << std::endl;
#endif
            for ( int i = 0; i < buf_sp.size(); ++i )
              N_scat[static_cast<particle::species>(buf_sp[i])] = buf_data[i];
          }
        }
        for ( const auto& rdir : sf.toc_dir() ) {
          if ( rdir.find("rank") != 0 ) continue;
          sf.cd(rdir);
          int r = sf.read1<int>("r");
#ifdef PIC_DEBUG
          lgr::file << "LDCKPT rank = " << r << std::endl;
          lgr::file << silo::errmsg() << std::endl;
#endif
          for ( auto i : sps ) {
            auto sp = static_cast<particle::species>(i);
            p_ckpt.load( sf,  properties[sp], particles[sp], myrank + r, ens_size );
          }
          sf.cd("..");
        }
        sf.cd("..");
      }
    }

    return checkpoint_ts;
  }


  template < int DGrid, typename R, template < typename > class S >
  std::string save_tracing( std::string prefix, const int num_parts,
                            const std::optional<dye::Ensemble<DGrid>>& ens_opt,
                            int timestep,
                            const particle::map<particle::array<R,S>>& particles
                            ) {
    particle::array<R,S> ptcs_tr {};
    if ( ens_opt ) {
      for ( auto sp : particles ) {
        for ( const auto& ptc : particles[sp] ) {
          if ( !ptc.is(particle::flag::exist) || !ptc.is(particle::flag::traced) ) continue;
          ptcs_tr.push_back(ptc);
        }
      }
    }
    bool is_ignore = !ens_opt || ( ptcs_tr.size() == 0 );
    int key = 0;
    if ( !is_ignore ) key = ens_opt -> label();
    // use ensemble label to put processes from same ensemble together
    auto active = mpi::world.split( {is_ignore}, key );
    if ( is_ignore ) return "";
    const auto& ens = *ens_opt;

    char str_ts [10];
    sprintf(str_ts, "%06d\0", timestep);

    prefix = prefix + "/tracing/timestep" + str_ts;

    silo::Pmpio pmpio;

    int num_per_part = active->size() / num_parts + ( active->size() % num_parts != 0 );

    {
      pmpio.filename = prefix + "/part" + std::to_string( active->rank() / num_per_part ) + ".silo";
      pmpio.dirname = "/ensemble" + std::to_string(ens.label());

      pmpio.comm = active -> split( active->rank() / num_per_part );
      fs::mpido(*active, [&](){fs::create_directories(prefix);} );
    }

    active->barrier();

    pmpio([&]( auto& dbfile )
          {
            // write global data
            if ( !dbfile.exists("/timestep") ) {
              using T = std::underlying_type_t<particle::species>;
              std::vector<T> sps;
              for ( auto sp : particles )
                sps.push_back( static_cast<T>(sp) );
              dbfile.write( "/species", sps  );

              dbfile.write( "/timestep", timestep );
            }

            dbfile.mkcd( "rank" + std::to_string(ens.intra.rank()) );
            ParticleArrayCkpt<R, S> p_ckpt;
            { // write process specific data
              int r = ens.intra.rank();
              dbfile.write( "r", r );
              {
                const auto& ptcs = ptcs_tr;
                particle::load_t num = ptcs.size();
                // NOTE: explicitly write load 1) to signal that this species in principle can exist and 2) because in Silo array length has type int so inferring from it might in some extreme case is wrong and hard to debug
                dbfile.write("load", num);
                if ( num != 0 ) {
                  constexpr auto DPtc = S<R>::Dim;
                  for ( int i = 0; i < DPtc; ++i ) {
                    dbfile.write("q"+std::to_string(i+1), p_ckpt.q(ptcs, i).data(), {p_ckpt.q(ptcs, i).size()} );
                    dbfile.write("p"+std::to_string(i+1), p_ckpt.p(ptcs, i).data(), {p_ckpt.p(ptcs, i).size()} );
                  }
                  dbfile.write( "frac", p_ckpt.frac(ptcs).data(), {p_ckpt.frac(ptcs).size()} );
                  dbfile.write( "state", p_ckpt.state(ptcs).data(), {p_ckpt.state(ptcs).size()} );
                }
              }

            }
            dbfile.cd("..");
          }
      );

    return prefix;
  }

  // template < int DGrid, typename R, template < typename > class S >
  // int load_tracing( std::string dir, particle::array<R,S>& particles) {
  //   int timestep = 0;
  //   using T = std::underlying_type_t<particle::species>;

  //   // rank0 read data
  //   if ( mpi::world.rank() == 0 ) {
  //     auto dir_itr = fs::directory_iterator(dir);
  //     auto dbfile = silo::open(*dir_itr, silo::Mode::Read);
  //     assert( dbfile.var_exists("/timestep") );

  //     dbfile.read("/timestep", &timestep);
  //   }
  //   // rank0 broadcast data
  //   {
  //     int buf[1] = {timestep};
  //     mpi::world.broadcast( 0, buf, 1 );
  //     if ( mpi::world.rank() != 0 ) {
  //       timestep = buf[0];
  //     }
  //   }

  //   // ----------

  //   const int myrank = ens_opt->intra.rank();

  //   ParticleArrayCkpt<R, S> p_ckpt;
  //   for ( auto f : fs::directory_iterator(dir) ) {
  //     auto sf = silo::open(f, silo::Mode::Read);
  //     for ( const auto& dname : sf.toc_dir() ) {
  //       if ( dname.find("ensemble") != 0 ) continue;
  //       sf.cd(dname);
  //       for ( const auto& rdir : sf.toc_dir() ) {
  //         if ( rdir.find("rank") != 0 ) continue;
  //         sf.cd(rdir);
  //         int r = sf.read1<int>("r");

  //         {
  //           auto N = dbfile.read1<particle::load_t>("load");

  //           if ( 0 != N ) {
  //             // auto[ from, num ] = dye::scatter_load(N, receiver_idx, num_receivers);
  //             int from = 0;
  //             auto num = N;
  //             auto size_old = ptcs.size();
  //             ptcs.resize(size_old + num);

  //             constexpr auto DPtc = S<T>::Dim;
  //             apt::array<int, 3> slice { from, num, 1 };
  //             {
  //               for ( int i = 0; i < DPtc; ++i ) {
  //                 dbfile.readslice("q"+std::to_string(i+1), {slice}, particles._q[i].data() + size_old );
  //                 dbfile.readslice("p"+std::to_string(i+1), {slice}, particles._p[i].data() + size_old );
  //               }
  //               dbfile.readslice("frac", {slice}, particles._frac.data() + size_old );
  //               dbfile.readslice("state", {slice}, particles._state.data() + size_old );
  //             }
  //           }
  //         }

  //         sf.cd("..");
  //       }
  //       sf.cd("..");
  //     }
  //   }

  //   return timestep;
  // }
}
