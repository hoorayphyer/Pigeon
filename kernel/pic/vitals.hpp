#ifndef _PIC_VITALS_HPP_
#define _PIC_VITALS_HPP_

#include "particle/load_type.hpp"
#include "particle/properties.hpp"
#include "apt/csi.hpp"
#include "timer/timer.hpp"
#include <fstream>

namespace pic {
  template < int DGrid,
             typename T,
             template < typename > class S
             >
  void check_vitals( std::string filename, int timestep, const dye::Ensemble<DGrid>& ens,
                     const std::optional<mpi::CartComm>& cart,
                     const particle::map<particle::array<T,S>>& particles,
                     const particle::map<particle::load_t>& N_scat ) {
    using namespace particle;
    static int counter = 0;
    constexpr int interval = 40;
    static tmr::Timestamp stopwatch;

    std::vector<load_t> buffer;
    {
      for ( auto sp : particles ) buffer.push_back(particles[sp].size());
      for ( auto sp : N_scat ) buffer.push_back(N_scat[sp]);
      ens.reduce_to_chief( mpi::by::SUM, buffer.data(), buffer.size() );
    }
    if ( !cart ) return;

    buffer.push_back(ens.size());
    cart->template reduce<true>( mpi::by::SUM, 0, buffer.data(), buffer.size() );

    if ( cart->rank() != 0 ) return;

    auto* p1 = buffer.data();
    auto* p2 = buffer.data() + particles.size();

    {
      std::ofstream out(filename, std::ios_base::app);
      if ( counter % interval == 0 ) {
        out << "timestep|\tnprocs|\tlapse/hr|\tTotal load|\tcumulative new ptcs from scattering" << std::endl;
        out << "species ordering : ";
        for ( auto sp : N_scat )
          out << properties[sp].nickname << " ";
        out << std::endl;
        counter = 0;
      }
      ++counter;
      float lps = stopwatch.lapse().in_units_of("s").val() / 3600.0;
      out << timestep << "|\t" << buffer.back() << "|\t" << lps <<  "|\t";
      for ( int i = 0; i < particles.size(); ++ i )
        out << apt::csi(p1[i]) << " ";
      out << "|\t";
      for ( int i = 0; i < N_scat.size(); ++i )
        out << apt::csi(p2[i]) << " ";
      out << std::endl;
      out.close();
    }

    return;
  }
}

#endif
