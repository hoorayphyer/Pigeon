#include "ckpt/checkpoint.hpp"
#include "mpipp/mpi++.hpp"
#include "silopp/silo++.hpp"
#include "silopp/pmpio.hpp"
#include "filesys/filesys.hpp"
#include "particle/properties.hpp"
#include "particle/load_type.hpp"
#include "field/field.hpp"
#include <cassert>
#include "dye/dynamic_balance.hpp"
#include "dye/scatter_load.hpp"

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
        dbfile.read( n, f._comps[comp].data() );
      }
    }
  };

  template < typename T, template < typename > class PtcSpecs >
  struct ParticleArrayCkpt {
    void save( silo::file_t& dbfile, const particle::species& sp, const particle::array<T,PtcSpecs>& ptcs ) {
      dbfile.mkcd(particle::properties.at(sp).name);
      particle::load_t num = ptcs.size();
      // NOTE: explicitly write load 1) to signal that this species in principle can exist and 2) because in Silo array length has type int so inferring from it might in some extreme case is wrong and hard to debug
      dbfile.write("load", num);
      if ( num != 0 ) {
        constexpr auto DPtc = PtcSpecs<T>::Dim;
        for ( int i = 0; i < DPtc; ++i ) {
          dbfile.write("q"+std::to_string(i+1), ptcs._q[i].data(), {ptcs._q[i].size()} );
          dbfile.write("p"+std::to_string(i+1), ptcs._p[i].data(), {ptcs._p[i].size()} );
        }
        dbfile.write( "state", ptcs._state.data(), {ptcs._state.size()} );
      }
      dbfile.cd("..");
    }

    void load( silo::file_t& dbfile, const particle::species& sp, particle::array<T,PtcSpecs>& ptcs, int receiver_idx, int num_receivers ) {
      dbfile.cd( particle::properties.at(sp).name );
      auto N = dbfile.read1<particle::load_t>("load");

      if ( 0 != N ) {
        auto[ from, num ] = dye::scatter_load(N, receiver_idx, num_receivers);
        auto size_old = ptcs.size();
        ptcs.resize(size_old + num);

        constexpr auto DPtc = PtcSpecs<T>::Dim;
        apt::array<int, 3> slice { from, num, 1 };
        {
          for ( int i = 0; i < DPtc; ++i ) {
            dbfile.readslice("q"+std::to_string(i+1), {slice}, ptcs._q[i].data() + size_old );
            dbfile.readslice("p"+std::to_string(i+1), {slice}, ptcs._p[i].data() + size_old );
          }
          dbfile.readslice("state", {slice}, ptcs._state.data() + size_old );
        }
      }
      dbfile.cd("..");
    }

  };
}

namespace ckpt {
  // TODO save random seed
  template < int DGrid,
             typename Real,
             template < typename > class PtcSpecs >
  std::string save_checkpoint( std::string prefix, const int num_parts,
                               const std::optional<dye::Ensemble<DGrid>>& ens_opt,
                               int timestep,
                               const field::Field<Real, 3, DGrid>& E,
                               const field::Field<Real, 3, DGrid>& B,
                               const particle::map<particle::array<Real,PtcSpecs>>& particles
                               ) {
    bool is_idle = !ens_opt;
    int key = 0;
    if ( !is_idle ) key = ens_opt -> label();
    // use ensemble label to put same ensemble processes together
    auto active = mpi::world.split( {is_idle}, key );
    if ( is_idle ) return "";
    const auto& ens = *ens_opt;

    char str_ts [10];
    sprintf(str_ts, "%06d\0", timestep);

    prefix = prefix + "/checkpoints/checkpoint_timestep" + str_ts;

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
            { // write global data
              if ( !dbfile.exists("/timestep") ) {
                dbfile.write( "/timestep", timestep );
              }
              if ( !dbfile.exists("/cartesian_partition") ) {
                dbfile.write( "/cartesian_partition", ens.cart_dims.begin(), {DGrid} );
              }
              if ( !dbfile.exists("/species") ) {
                using T = std::underlying_type_t<particle::species>;
                std::vector<T> sps;
                for ( const auto& [sp, ptcs] : particles )
                  sps.push_back( static_cast<T>(sp) );
                dbfile.write( "/species", sps  );
              }
            }

            { // write ensemble-wise data
              if ( !dbfile.exists("label") ) {
                dbfile.write( "label", ens.label() );
              }

              if ( !dbfile.exists("cartesian_coordinates") ) {
                dbfile.write( "cartesian_coordinates", ens.cart_coords.begin(), {DGrid} );
              }

              std::vector<particle::load_t> loads;
              loads.reserve(particles.size());
              for ( const auto& [sp, ptcs] : particles )
                loads.push_back(ptcs.size());

              ens.reduce_to_chief( mpi::by::SUM, loads.data(), loads.size() );

              FieldCkpt<Real,DGrid> ckpt;
              if ( ens.intra.rank() == ens.chief ) {
                int idx = 0;
                for ( const auto& [sp, ptcs] : particles )
                  dbfile.write( particle::properties[sp].name + "_load", loads[idx++] );

                ckpt.save(dbfile, "E", E);
                ckpt.save(dbfile, "B", B);
              }
            }

            dbfile.mkcd( "rank" + std::to_string(ens.intra.rank()) );
            { // write process specific data
              int r = ens.intra.rank();
              dbfile.write( "r", r );
              ParticleArrayCkpt<Real, PtcSpecs> ckpt;
              for ( const auto& [sp, ptcs] : particles ) {
                ckpt.save( dbfile, sp, ptcs );
              }
            }
            dbfile.cd("..");
          }
      );

    return prefix;

  }

  template < int DGrid,
             typename Real,
             template < typename > class PtcSpecs >
  int load_checkpoint( std::string dir,
                       std::optional<dye::Ensemble<DGrid>>& ens_opt,
                       const std::optional<mpi::CartComm>& cart_opt,
                       field::Field<Real, 3, DGrid>& E,
                       field::Field<Real, 3, DGrid>& B,
                       particle::map<particle::array<Real,PtcSpecs>>& particles,
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
      assert( dbfile.var_exists("/cartesian_partition") );
      assert( dbfile.var_exists("/species") );

      dbfile.read("/timestep", &checkpoint_ts);
      {
        int ndims = dbfile.var_length("/cartesian_partition");
        std::vector<int> parts(ndims);
        dbfile.read("/cartesian_partition", parts.data());
        num_ens = 1;
        for ( auto x : parts ) num_ens *= x;
      }
      {
        num_sps = dbfile.var_length("/species");
        sps.resize( num_sps );
        dbfile.read("/species", sps.data() );
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
                std::string entry = particle::properties[sp].name + "_load";
                assert( sf.var_exists(entry) );
                myload += sf.read1<particle::load_t>(entry);
              }
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
    {
      const int mylabel = ens_opt->label();
      const int myrank = ens_opt->intra.rank();
      const auto ens_size = ens_opt->intra.size();

      FieldCkpt<Real,DGrid> f_ckpt;
      ParticleArrayCkpt<Real, PtcSpecs> p_ckpt;
      for ( auto f : fs::directory_iterator(dir) ) {
        auto sf = silo::open(f, silo::Mode::Read);
        for ( const auto& dname : sf.toc_dir() ) {
          if ( dname.find("ensemble") != 0 ) continue;
          int l = sf.read1<int>(dname+"/label");
          if ( l != mylabel ) continue;
          sf.cd(dname);
          if ( 0 == myrank ) {
            f_ckpt.load( sf, "E", E );
            f_ckpt.load( sf, "B", B );
          }
          for ( const auto& rdir : sf.toc_dir() ) {
            if ( rdir.find("rank") != 0 ) continue;
            sf.cd(rdir);
            int r = sf.read1<int>("r");
            for ( auto i : sps ) {
              auto sp = static_cast<particle::species>(i);
              p_ckpt.load( sf, sp, particles.at(sp), myrank + r, ens_size );
            }
            sf.cd("..");
          }
          sf.cd("..");
        }
      }
    }

    return checkpoint_ts;
  }
}
