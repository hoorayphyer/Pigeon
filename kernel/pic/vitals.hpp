#ifndef _PIC_VITALS_HPP_
#define _PIC_VITALS_HPP_

#include "particle/properties.hpp"
#include "timer/timer.hpp"
#include <fstream>

namespace pic {
  template < int DGrid,
             typename T,
             template < typename > class S
             >
  void check_vitals( std::string filename, T t_phys, const dye::Ensemble<DGrid>& ens,
                     const std::optional<mpi::CartComm>& cart,
                     const particle::map<particle::Properties>& properties,
                     const particle::map<particle::array<T,S>>& particles,
                     const particle::map<T>& N_scat ) {
    using namespace particle;
    static int counter = 0;
    constexpr int interval = 40; // how often to print the table description row
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

    // significant only on world rank 0
    static double t_phys_prev = -1e-6;
    static std::vector<double> num_ptcs_prev (particles.size(), 0);
    static std::vector<double> num_scat_prev (particles.size(), 0);

    auto* p1 = buffer.data();
    auto* p2 = buffer.data() + particles.size();

    {
      std::ofstream out(filename, std::ios_base::app);
      if ( counter % interval == 0 ) {
        out << "t_phys|\tnprocs|\tlapse/hr|\tTotal load(rate)|\tscattering creation rate" << std::endl;
        out << "species ordering : ";
        for ( auto sp : N_scat )
          out << properties[sp].nickname << " ";
        out << std::endl;
        counter = 0;
      }
      ++counter;
      out << apt::fmt("%.2f",t_phys) << "|\t" << buffer.back() << "|\t" << apt::fmt("%8.2f", stopwatch.lapse().in_units_of("s").val() / 3600.0) <<  "|\t";
      for ( int i = 0; i < particles.size(); ++ i ) {
        out << apt::fmt("%.2e", p1[i]) << "(" << apt::fmt( "%.2e", (p1[i] - num_ptcs_prev[i]) / (t_phys- t_phys_prev) ) << ") ";
        num_ptcs_prev[i] = p1[i];
      }
      out << "|\t";
      for ( int i = 0; i < N_scat.size(); ++i ) {
        out <<  apt::fmt( "%.2e", (p2[i] - num_scat_prev[i]) / (t_phys- t_phys_prev) ) << " ";
        num_scat_prev[i] = p2[i];
      }
      out << std::endl;
      out.close();

      t_phys_prev = t_phys;
    }

    return;
  }
}

#endif
