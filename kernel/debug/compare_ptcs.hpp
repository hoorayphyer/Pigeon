#ifndef _DEBUG_COMPARE_PTCS_HPP_
#define _DEBUG_COMPARE_PTCS_HPP_

#include "particle/sorter.hpp"
#include "mpipp/mpi++.hpp"

namespace debug {
  // compare two sets of particles which may or may not be distributed over processes modulo empty particles.
  template < typename T, template < typename > class S >
  std::array<std::size_t, 5> compare_particles( particle::array<T, S>&& ptcs1, particle::array<T, S>&& ptcs2, const mpi::Comm& comm ) {
    using namespace particle;
    auto is_same =
      []( const auto& ptc1, const auto& ptc2 ) {
        for ( int i = 0; i < S<T>::Dim; ++i ) {
          if ( ptc1.q(i) != ptc2.q(i) ) return false;
          if ( ptc1.p(i) != ptc2.p(i) ) return false;
        }
        if ( ptc1.state() != ptc2.state() ) return false;
        return true;
      };

    auto count_ptcs =
      [] ( const auto& ptcs ) {
        std::size_t num = 0;
        for ( const auto& x : ptcs ) num += x.is(flag::exist);
        return num;
      };

    const std::size_t num_ptcs1_before = count_ptcs(ptcs1);
    const std::size_t num_ptcs2_before = count_ptcs(ptcs2);

    const int r = comm.rank();
    const int size = comm.size();

    std::size_t count_matched = 0;

    int n = comm.size();

    do {
      for ( auto x : ptcs1 ) { // TODOL semantics
        if ( !x.is(flag::exist) ) continue;
        for ( auto y : ptcs2 ) {
          if ( !y.is(flag::exist) ) continue;
          if ( !is_same(x,y) ) continue;
          ++count_matched;
          x.reset( flag::exist );
          y.reset( flag::exist );
        }
      }
      sort( ptcs1 );
      sort( ptcs2 );
      // pass ptcs1 to next process for further cancellation
      {
        std::size_t count_s =  ptcs1.size();
        std::vector<Particle<T,S>> buf_s ( count_s );
        std::vector<Particle<T,S>> buf_r;
        for ( int i = 0; i < count_s; ++i )
          buf_s[i] = std::move(ptcs1[i]);

        std::vector<mpi::Request> reqs(2);
        reqs[0] = comm.Isend( (r + 1) % size, 147, buf_s.data(), buf_s.size() );

        std::size_t count_r = comm.probe( (r + size - 1) % size, 147, buf_r.data() );
        buf_r.resize(count_r);
        reqs[1] = comm.Irecv( (r + size - 1) % size, 147, buf_r.data(), count_r );
        mpi::waitall(reqs);

        for ( int i = 0; i < count_r; ++i )
          ptcs1[i] = std::move(buf_r[i]);
      }

    } while ( --n );

    return { num_ptcs1_before, num_ptcs2_before, count_matched, count_ptcs(ptcs1), count_ptcs(ptcs2)};

  }
}

#endif
