#include "dye/scatter_load.hpp"
#include "particle/array.hpp"
#include "particle/load_type.hpp"
#include "silopp/silo++.hpp"

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
#if PIC_DEBUG
        lgr::file << "LDCKPT Field " << n << ", length of ckpt data = " << dbfile.var_length(n) << ", length of container = " << f._comps[comp].size() << std::endl;
#endif
        dbfile.read( n, f._comps[comp].data() );
#if PIC_DEBUG
        lgr::file << silo::errmsg() << std::endl;
#endif
      }
    }
  };

  template < typename T, template < typename > class S >
  struct ParticleArrayCkpt {
    void save( silo::file_t& dbfile, const std::string& varname, const particle::array<T,S>& ptcs ) {
      dbfile.mkcd(varname);
      particle::load_t num = ptcs.size();
      // NOTE: explicitly write load 1) to signal that this species in principle can exist and 2) because in Silo array length has type int so inferring from it might in some extreme case is wrong and hard to debug
      dbfile.write("load", num);
      if ( num != 0 ) {
        constexpr auto DPtc = S<T>::Dim;
        for ( int i = 0; i < DPtc; ++i ) {
          dbfile.write("q"+std::to_string(i+1), ptcs.qs(i).data(), {ptcs.qs(i).size()} );
          dbfile.write("p"+std::to_string(i+1), ptcs.ps(i).data(), {ptcs.ps(i).size()} );
        }
        dbfile.write( "frac", ptcs.fracs().data(), {ptcs.fracs().size()} );
        dbfile.write( "state", ptcs.states().data(), {ptcs.states().size()} );
      }
      dbfile.cd("..");
    }

    void load( silo::file_t& dbfile, const std::string& varname, particle::array<T,S>& ptcs, int receiver_idx=0, int num_receivers=1 ) {
      dbfile.cd( varname );
      auto N = dbfile.read1<particle::load_t>("load");
#if PIC_DEBUG
      lgr::file << "LDCKPT SPECIES " << varname << ", Load = " << N << std::endl;
#endif

      if ( 0 != N ) {
        auto[ from, num ] = dye::scatter_load(N, receiver_idx, num_receivers);
        auto size_old = ptcs.size();
        ptcs.resize(size_old + num);

        constexpr auto DPtc = S<T>::Dim;
        silo::Slice slice { from, from+num, 1 };
#if PIC_DEBUG
        lgr::file << "LDCKPT Slice from = " << from << ", length = " << num << ", stride = " << 1 << std::endl;
#endif
        {
          for ( int i = 0; i < DPtc; ++i ) {
            dbfile.readslice("q"+std::to_string(i+1), {slice}, ptcs.qs(i).data() + size_old );
#if PIC_DEBUG
            lgr::file << silo::errmsg() << std::endl;
#endif
            dbfile.readslice("p"+std::to_string(i+1), {slice}, ptcs.ps(i).data() + size_old );
#if PIC_DEBUG
            lgr::file << silo::errmsg() << std::endl;
#endif
          }
          dbfile.readslice("frac", {slice}, ptcs.fracs().data() + size_old );
          dbfile.readslice("state", {slice}, ptcs.states().data() + size_old );
#if PIC_DEBUG
          lgr::file << silo::errmsg() << std::endl;
#endif
        }
      }
      dbfile.cd("..");
    }

  };
}
