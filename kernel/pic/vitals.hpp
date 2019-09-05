#ifndef _PIC_VITALS_HPP_
#define _PIC_VITALS_HPP_

#include "particle/properties.hpp"
#include "timer/timer.hpp"
#include "msh/mesh_shape_interplay_impl.hpp" // WeightFinder
#include <fstream>
#include <cstdio>

namespace pic {
  std::string sci( double x ) {
    char str[20];
    std::snprintf(str, 20, "%.5e", x);
    return {str};
  }

  template < int DGrid,
             typename T,
             template < typename > class S
             >
  void check_vitals( std::string filename, int timestep, const dye::Ensemble<DGrid>& ens,
                     const std::optional<mpi::CartComm>& cart,
                     const particle::map<particle::array<T,S>>& particles,
                     const particle::map<double>& N_scat ) {
    using namespace particle;
    static int counter = 0;
    constexpr int interval = 40;
    static tmr::Timestamp stopwatch;

    std::vector<double> buffer;
    {
      for ( auto sp : particles ) {
        double num = 0.0;
        for ( const auto& ptc : particles[sp] ) {
          if ( ptc.is(flag::exist) ) num += ptc.frac();
        }
        buffer.push_back(num);
      }
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
        out << sci(p1[i]) << " ";
      out << "|\t";
      for ( int i = 0; i < N_scat.size(); ++i )
        out << sci(p2[i]) << " ";
      out << std::endl;
      out.close();
    }

    return;
  }
}

#endif
