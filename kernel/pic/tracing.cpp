#include "pic/tracing.hpp"
#include "pic.hpp"

#include "dye/ensemble.hpp"
#include "filesys/filesys.hpp"
#include "mpipp/mpi++.hpp"
#include "particle/properties.hpp"
#include "particle/load_type.hpp"
#include "silopp/pmpio.hpp"
#include "silopp/silo++.hpp"

namespace particle {
  using R = pic::real_t;

  template <typename T>
  using S = particle::Specs<T>;

  constexpr int DGrid = pic::DGrid;

  template <>
  std::unordered_map<species, std::size_t>
  TracingManager::parse_tracing(const map<array<R, S>> &particles) {
    std::unordered_map<species, std::size_t> res;
    const auto &comm = mpi::world;

    for ( auto sp : particles ) {
      std::vector<std::size_t> counts(comm.size(),{});

      for ( const auto& p : particles[sp] ) {
        std::size_t id = p.template get<pid>();
        std::size_t sn = apt::getbits<1 + Nbits_wr, Nbits_sn>(id);
        std::size_t wr = apt::getbits<1, Nbits_wr>(id);
        if ( wr < comm.size() )
          counts[wr] = std::max(counts[wr], sn);
      }

      comm.allreduce<mpi::IN_PLACE>(mpi::by::MAX, counts.data(), counts.size() );
      res[sp] = counts[comm.rank()];
    }

    return res;
  }


  template <>
  std::string
  TracingManager::save_tracing(std::string prefix, const int num_parts,
                               const std::optional<dye::Ensemble<DGrid>> &ens_opt,
                               int timestep,
                               const map<array<R, S>> &particles) {
    particle::array<R,S> ptcs_tr {};
    if ( ens_opt ) {
      for ( auto sp : particles ) {
        for ( const auto& ptc : particles[sp] ) {
          if ( !ptc.is(flag::exist) || !is_traced(ptc) ) continue;
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
                    dbfile.write("q"+std::to_string(i+1), ptcs.qs(i).data(), {ptcs.qs(i).size()} );
                    dbfile.write("p" + std::to_string(i + 1), ptcs.ps(i).data(), {ptcs.ps(i).size()});
                  }
                  dbfile.write( "frac", ptcs.fracs().data(), {ptcs.fracs().size()} );
                  dbfile.write( "state", ptcs.states().data(), {ptcs.states().size()} );
                }
              }

            }
            dbfile.cd("..");
          }
      );

    return prefix;
  }
}
