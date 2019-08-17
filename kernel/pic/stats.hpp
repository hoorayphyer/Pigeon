#ifndef _PARTICLE_STATS_HPP_
#define _PARTICLE_STATS_HPP_

#include "particle/load_type.hpp"
#include "particle/properties.hpp"
#include <fstream>
#include <limits>

namespace particle {
  template < int DGrid,
             typename T,
             template < typename > class S
             >
  void statistics( std::string filename, int timestep, const dye::Ensemble<DGrid>& ens, const std::optional<mpi::CartComm>& cart, const map<array<T,S>>& particles) {
    static bool first_time = true;
    std::vector<load_t> buffer; // layout of buffer: ens_size, sum_sp1, max_sp1, min_sp1, sum_sp2, max_sp2, min_sp2...
    {
      if ( ens.is_chief() ) buffer.push_back(ens.intra.size());
      for ( const auto& [sp, ptcs] : particles ) {
        load_t load = ptcs.size();
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
        for ( const auto& [sp, ignore ] : particles ) {
          ave[sp] = 0;
          max[sp] = std::numeric_limits<load_t>::min();
          min[sp] = std::numeric_limits<load_t>::max();
          sp2idx[sp] = idx++;
        }
      }
      for( int i = 0; i < cart->size(); i += stride ) {
        nprocs += stats[i * stride];
        for ( const auto& [sp, ignore ] : particles ) {
          ave[sp] += stats[ i * stride + 1 + 3 * sp2idx[sp] ];
          max[sp] = std::max<load_t>( max[sp], stats[ i * stride + 1 + 3 * sp2idx[sp] + 1 ] );
          min[sp] = std::min<load_t>( min[sp], stats[ i * stride + 1 + 3 * sp2idx[sp] + 2 ] );
        }
      }
      for ( const auto& [sp, ignore ] : particles ) ave[sp] /= nprocs;
    }

    {
      std::ofstream out(filename, std::ios_base::app);
      if ( first_time ) {
        out << "timestep|\tnprocs";
        for ( const auto& [ sp, ignore ] : particles )
          out << "|\t" << properties[sp].name;
        out << std::endl;
        out << "particle field : ave, max, min" << std::endl;
        first_time = false;
      }
      out << timestep << "|\t" << nprocs;
      for ( const auto& [ sp, ignore ] : particles )
        out << "|\t" << ave[sp] << ", " << max[sp] << ", " << min[sp];
      out << std::endl;
      out.close();
    }

    return;
  }
}

#endif
