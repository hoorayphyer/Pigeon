#include "ckpt/checkpoint.hpp"
#include "mpipp/mpi++.hpp"
#include "silopp/silo++.hpp"
#include "silopp/pmpio.hpp"
#include "filesys/filesys.hpp"
#include "particle/properties.hpp"

// auxilliary
namespace ckpt {
  template < typename T, int DGrid >
  struct FieldCkpt {
    template < int DField >
    void save( silo::file_t& dbfile, const std::string& name, const field::Field<T,DField,DGrid>& f ) {
      std::vector<int> dims(DGrid);
      for ( int comp = 0; comp < DField; ++comp ) {
        auto n = name + ( DField == 1 ? "" : std::to_string(comp+1) );
        dbfile.write( n, f._comps[comp] );
      }
    }

    template < int DField >
    void load( silo::file_t& dbfile, const std::string& name, const field::Field<T,DField,DGrid>& f ) {}
  };

  template < typename T, template < typename > class PtcSpecs >
  struct ParticleArrayCkpt {
    void save( silo::file_t& dbfile, const particle::species& sp, const particle::array<T,PtcSpecs>& ptcs ) {
      if ( ptcs.size() == 0 ) return;
      constexpr auto DPtc = PtcSpecs<T>::Dim;
      dbfile.mkcd(particle::properties.at(sp).name);
      {
        for ( int i = 0; i < DPtc; ++i ) {
          dbfile.write("q"+std::to_string(i+1), ptcs._q[i].data(), {ptcs._q[i].size()} );
          dbfile.write("p"+std::to_string(i+1), ptcs._p[i].data(), {ptcs._p[i].size()} );
        }
        dbfile.write( "state", ptcs._state.data(), {ptcs._state.size()} );
      }
      dbfile.cd("..");
    }

    void load( silo::file_t& dbfile, const particle::species& sp, const particle::array<T,PtcSpecs>& ptcs ) {
    }

  };
}

namespace ckpt {
  template < int DGrid,
             typename Real,
             template < typename > class PtcSpecs >
  void save_checkpoint( std::string prefix, const int num_parts,
                      const std::optional<dye::Ensemble<DGrid>>& ens_opt,
                      int timestep,
                      const field::Field<Real, 3, DGrid>& E,
                      const field::Field<Real, 3, DGrid>& B,
                      const particle::map<particle::array<Real,PtcSpecs>>& particles
                      ) {
    bool is_idle = !ens_opt;
    auto active = mpi::world.split( {is_idle} );
    if ( is_idle ) return;
    const auto& ens = *ens_opt;

    char str_ts [10];
    sprintf(str_ts, "%06d\0", timestep);

    prefix = prefix + "/checkpoints/checkpoint_timestep" + str_ts;

    silo::Pmpio pmpio;
    {
      pmpio.filename = prefix + "/part" + std::to_string( active->rank() % num_parts ) + ".silo";
      pmpio.dirname = "/ensemble" + std::to_string(ens.label());

      pmpio.comm = active -> split( active->rank() % num_parts );
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
            }

            { // write ensemble-wise data
              if ( !dbfile.exists("cartesian_coordinates") ) {
                dbfile.write( "cartesian_coordinates", ens.cart_coords.begin(), {DGrid} );
              }
              FieldCkpt<Real,DGrid> ckpt;
              if ( ens.intra.rank() == ens.chief ) {
                ckpt.save(dbfile, "E", E);
                ckpt.save(dbfile, "B", B);
              }
            }

            dbfile.mkcd( "rank" + std::to_string(ens.intra.rank()) );
            { // write process specific data
              ParticleArrayCkpt<Real, PtcSpecs> ckpt;
              for ( const auto& [sp, ptcs] : particles ) {
                ckpt.save( dbfile, sp, ptcs );
              }
            }
            dbfile.cd("..");
          }
      );


  }
}
