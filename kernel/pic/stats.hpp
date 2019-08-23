#ifndef _PARTICLE_STATS_HPP_
#define _PARTICLE_STATS_HPP_

#include "particle/load_type.hpp"
#include "particle/properties.hpp"
#include "apt/csi.hpp"
#include <fstream>
#include <limits>

namespace particle {
  template < int DGrid,
             typename T,
             template < typename > class S
             >
  void statistics( std::string filename, int timestep, const dye::Ensemble<DGrid>& ens, const std::optional<mpi::CartComm>& cart, const map<array<T,S>>& particles, map<load_t>& N_scat ) {
    static bool first_time = true;
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
          ave[sp] = 0;
          max[sp] = std::numeric_limits<load_t>::min();
          min[sp] = std::numeric_limits<load_t>::max();
          sp2idx[sp] = idx++;
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
      if ( first_time ) {
        out << "timestep|\tnprocs";
        for ( auto sp : particles )
          out << "|\t" << properties[sp].name;
        out << "|\t" << "cumulative new ptcs from scattering" << std::endl;
        out << "particle field : ave, max, min" << std::endl;
        out << "scattering : ";
        for ( auto sp : N_scat )
          out << properties[sp].nickname << " ";
        out << std::endl;
        first_time = false;
      }
      out << timestep << "|\t" << nprocs;
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
