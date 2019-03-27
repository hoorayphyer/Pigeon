#include "particle/migration.hpp"
#include "apt/type_traits.hpp"
#include "parallel/mpi++.hpp"
#include <memory>

namespace particle :: impl {

  // NOTE Assume there is no empty particles.
  template < typename Buffer, typename F_LCR >
  auto lcr_sort( Buffer& buffer, const F_LCR& lcr ) noexcept {
    constexpr int L = 0;
    constexpr int C = 1;
    constexpr int R = 2;

    // ec points to end of ccc.., el points to end of lll..., and er points to the reverse end of ...rrr.
    // the upshot: ccc.. will be [0, ic), lll... will be [ic, ilr+1), rrr... will be [ilr+1, buffer.size()), edge case safe
    int ec = 0, el = 0, er = buffer.size() - 1;
    while ( el <= er ) {
      switch ( lcr( buffer[el] ) ) {
      case L : el++; break;
      case C :
        if ( ec != el ) std::swap( buffer[el], buffer[ec] );
        ec++; el++; break;
      case R :
        while ( el <= er && lcr( buffer[er] ) != R ) er--;
        if ( er != el ) std::swap( buffer[el], buffer[er] );
        er--;
        break;
      }
    }

    return apt::array<const int, 3>{ ec, el, buffer.size() };
  }


  template < typename Ptc, typename F_LCR >
  void migrate_1dim ( std::vector<Ptc>& buffer,
                      const apt::pair<std::optional<mpi::InterComm>>& intercomms,
                      const F_LCR& lcr, unsigned int shift ) {
    // sort order is center | left | right | empty. Returned are the delimiters between these catogories
    const auto begs = lcr_sort( buffer, lcr ); // begs = { begL, begR, begE_original };

    // NOTE buffer may be relocated in response to size growing. DO NOT store its pointers or references
    int begE_run = begs[2]; // running begin of empty particles in buffer
    for ( int lr = 0; lr < 2; ++lr ) {
      std::vector<mpi::Request> reqs;
      // sending
      if ( intercomms[lr] ) {
        const auto& send_comm = *intercomms[lr];
        int local_rank = send_comm.rank();
        int remote_dest = ( local_rank + shift ) % send_comm.remote_size();
        reqs.push_back( send_comm.Isend( remote_dest, 147, buffer.data() + begs[lr], begs[lr+1] - begs[lr] ) );
      }

      // receiving
      std::unique_ptr<Ptc[]> p_tmp(nullptr);
      int tot_num_recv = 0;
      if ( intercomms[1-lr] ) {
        const auto& recv_comm = *intercomms[1-lr];
        int local_rank = recv_comm.rank();
        int local_size = recv_comm.size();
        int remote_size = recv_comm.remote_size();
        std::vector<int> remote_srcs;
        std::vector<int> scan_recv_counts = {0}; // exclusive scan
        int src_rank = ( local_rank + local_size - (shift % local_size) ) % local_size;
        while ( src_rank < remote_size ) {
          remote_srcs.push_back(src_rank);
          scan_recv_counts.push_back( scan_recv_counts.back() + recv_comm.probe( src_rank, 147, buffer[0] ) );
          src_rank += local_size;
        }

        tot_num_recv = scan_recv_counts.back();
        // If recved more than space allows, store them in a temporary buffer then later merge with the primary buffer
        Ptc* p_recv = buffer.data() + begE_run;
        if ( tot_num_recv > buffer.capacity() - begE_run ) {
          p_tmp.reset( new Ptc [tot_num_recv] );
          p_recv = p_tmp.get();
        }

        for ( int i = 0; i < remote_srcs.size(); ++i ) {
          reqs.push_back( recv_comm.Irecv( remote_srcs[i], 147, p_recv + scan_recv_counts[i], scan_recv_counts[i+1] - scan_recv_counts[i] ) );
        }

      }

      mpi::waitall(reqs);

      // merge into buffer if needed.
      if ( intercomms[1-lr] && ( tot_num_recv > buffer.capacity() - begE_run ) ) {
        buffer.resize( begE_run + tot_num_recv );
        for ( int i = 0; i < tot_num_recv; ++i )
          buffer[begE_run + i] = p_tmp[i];
        p_tmp.reset(nullptr);
      }
      begE_run += tot_num_recv;
    }

    // erase sent particles by shifting
    for ( int i = 0; i < begE_run - begs[2]; ++i )
      buffer[begs[0] + i] = buffer[begs[2] + i];

    // NOTE It is essential to not have empty particles within buffer.size otherwise lcr_sort will fail.
    buffer.resize(begs[0] + begE_run - begs[2]);
  }
}

namespace particle {
  template < typename Vec, int DGrid, typename T >
  bool is_migrate( const apt::VecExpression<Vec,T>& q,
                   const apt::array< apt::pair<T>, DGrid>& borders ) noexcept {
    bool res = false;
    apt::foreach<0,DGrid>
      ( [&res]( const auto& x, const auto& bd ) noexcept {
          res |= ( x < bd[0] || x > bd[1] );
        }, q, borders );

    return res;
  }

  namespace impl {
    template < int I, typename T, int DPtc, typename state_t, int DGrid >
    inline void migrate( std::vector<cParticle<T,DPtc,state_t>>& buffer,
                         const apt::array< apt::pair<std::optional<mpi::InterComm>>, DGrid >& intercomms,
                         const apt::array< apt::pair<T>, DGrid >& borders,
                         unsigned int shift ) {
      auto&& lcr =
        [&bd=borders[I]]( const auto& ptc ) noexcept {
          return ( ptc.q()[I] >= bd[0] ) + ( ptc.q()[I] > bd[1] );
        };
      impl::migrate_1dim( buffer, intercomms[I], std::move(lcr), shift );
      if constexpr ( I > 0 )
                     migrate<I-1>( buffer, intercomms, borders, shift );
    }
  }

  template < typename T, int DPtc, typename state_t, int DGrid >
  void migrate ( std::vector<cParticle<T,DPtc,state_t>>& buffer,
                 const apt::array< apt::pair<std::optional<mpi::InterComm>>, DGrid >& intercomms,
                 const apt::array< apt::pair<T>, DGrid >& borders,
                 unsigned int pairing_shift ) {
    impl::migrate<DGrid-1>( buffer, intercomms, borders, pairing_shift );
  }

}
