#ifndef _PIC_VITALS_HPP_
#define _PIC_VITALS_HPP_

#include "particle/load_type.hpp"
#include "particle/properties.hpp"
#include "apt/csi.hpp"
#include "timer/timer.hpp"
#include <fstream>
#include <limits>

namespace pic {
  template < int DGrid,
             typename T,
             template < typename > class S
             >
  void check_vitals( std::string filename, int timestep, const dye::Ensemble<DGrid>& ens,
                     const std::optional<mpi::CartComm>& cart,
                     const particle::map<particle::array<T,S>>& particles,
                     particle::map<particle::load_t>& N_scat ) {
    using namespace particle;
    static int counter = 0;
    constexpr int interval = 40;
    static tmr::Timestamp stopwatch;

    std::vector<load_t> buffer;
    { // reduce N_scat to world rank 0
      ens.reduce_to_chief( mpi::by::SUM, N_scat.data().data(), N_scat.data().size() );
      if ( !ens.is_chief() ) {
        for( auto& x : N_scat.data() ) x = 0;
      }
      if ( cart ) {
        cart->template reduce<true>( mpi::by::SUM, 0, N_scat.data().data(), N_scat.data().size() );
        if ( cart->rank() != 0 ) {
          for( auto& x : N_scat.data() ) x = 0;
        }
      }
    }
    buffer.resize(0);
    { // layout of buffer: ens_size, sum_sp1, max_sp1, min_sp1, sum_sp2, max_sp2, min_sp2...
      if ( ens.is_chief() ) buffer.push_back(ens.intra.size());
      for ( auto sp : particles ) {
        load_t load = particles[sp].size();
        auto res = ens.intra.template reduce<false>( mpi::by::SUM, ens.chief, &load, 1 );
        if ( ens.is_chief() ) buffer.push_back((*res)[0]);
        res = ens.intra.template reduce<false>( mpi::by::MAX, ens.chief, &load, 1 );
        if ( ens.is_chief() ) buffer.push_back((*res)[0]);
        res = ens.intra.template reduce<false>( mpi::by::MIN, ens.chief, &load, 1 );
        if ( ens.is_chief() ) buffer.push_back((*res)[0]);
      }
      if ( !cart ) return;
    }

    auto stats_opt = cart->gather( 0, buffer.data(), buffer.size() );
    if ( cart->rank() != 0 ) return;

    const auto& stats = *stats_opt;
    auto stride = buffer.size();

    map<load_t> ave;
    map<load_t> max;
    map<load_t> min;
    int nprocs = 0;
    {
      map<int> sp2idx;
      {
        int idx = 0;
        for ( auto sp : particles ) {
          ave.insert(sp,0);
          max.insert(sp, std::numeric_limits<load_t>::min());
          min.insert(sp, std::numeric_limits<load_t>::max());
          sp2idx.insert(sp, idx++);
        }
      }
      for( int i = 0; i < cart->size(); ++i ) {
        nprocs += stats[i*stride];
        auto* s = &(stats[i*stride+1]);
        for ( auto sp : particles ) {
          ave[sp] += s[ 3 * sp2idx[sp] ];
          max[sp] = std::max<load_t>( max[sp], s[ 3 * sp2idx[sp] + 1 ] );
          min[sp] = std::min<load_t>( min[sp], s[ 3 * sp2idx[sp] + 2 ] );
        }
      }
      for ( auto sp : particles ) ave[sp] /= nprocs;
    }

    {
      std::ofstream out(filename, std::ios_base::app);
      if ( counter % interval == 0 ) {
        out << "timestep|\tnprocs|\tlapse/hr";
        for ( auto sp : particles )
          out << "|\t" << properties[sp].name;
        out << "|\t" << "cumulative new ptcs from scattering" << std::endl;
        out << "particle field : ave, max, min" << std::endl;
        out << "scattering : ";
        for ( auto sp : N_scat )
          out << properties[sp].nickname << " ";
        out << std::endl;
        counter = 0;
      }
      ++counter;
      float lps = stopwatch.lapse().in_units_of("s").val() / 3600.0;
      out << timestep << "|\t" << nprocs << "|\t" << lps;
      for ( auto sp : particles )
        out << "|\t" << apt::csi(ave[sp]) << ", " << apt::csi(max[sp]) << ", " << apt::csi(min[sp]);
      out << "|\t";
      for ( auto sp : N_scat )
        out << apt::csi(N_scat[sp]) << ", ";
      out << std::endl;
      out.close();
    }

    return;
  }
}

#endif
