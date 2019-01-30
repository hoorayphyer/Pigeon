#include "particle/migration.hpp"
#include "particle/particle_vector.hpp"
#include "parallel/mpi++.hpp"

namespace particle :: impl {

  // NOTE Assume there is no empty particles.
  template < typename T, std::size_t DPtc, typename F_LCR >
  auto lcr_sort( particle::vector<T,DPtc>& buffer, const F_LCR& lcr ) noexcept {
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
        if ( ec != el ) apt::swap( buffer[el], buffer[ec] );
        ec++; el++ break;
      case R :
        while ( el <= er && lcr( buffer[er] ) != R ) er--;
        if ( er != el ) apt::swap( buffer[el], buffer[er] );
        er--;
        break;
      }
    }

    return std::make_tuple( ec, el, buffer.size() );
  }


  template < typename T, std::size_t DPtc, typename F_LCR >
  void migrate_1dim ( particle::vector<T,DPtc>& buffer,
                      const std::array<int, 2>& neigh,
                      const F_LCR& lcr;
                      const mpi::Comm& comm ) {
    const auto& nL = std::get<0>(neigh);
    const auto& nR = std::get<1>(neigh);
    // sort order is center | left | right | empty. Returned are the delimiters between these catogories
    auto[begL, begR, begE] = lcr_sort( buffer, lcr );
    auto* ptr = buffer.data();
    int num_recv = 0;

    // TODO finish the comm code
    // send to left
    comm.send( ptr + begL, begR - begL, nL );
    num_recv += comm.recv( ptr + begE + num_recv, buffer.capacity() - begE - num_recv, nR );
    comm.wait();

    // send to right
    comm.send( ptr + begR, begE - begR, nR );
    num_recv += comm.recv( ptr + begE + num_recv, buffer.capacity() - begE - num_recv, nL );
    comm.wait();

    // erase sent particles by shifting
    for ( int i = 0; i < num_recv; ++i )
      ptr[begL + i] = ptr[begE + i];
    // erase outstanding sent particles
    buffer.erase( ptr + begL + num_recv, std::max(0, begE - begL - num_recv) ); // TODO erase( int from, std::size_t n );
    buffer.size = begL + num_recv;
  }
}

namespace particle {
  template < typename Tvt, std::size_t DGrid,
             typename Trl = apt::remove_cvref_t<Tvt>,
             typename T_q = Vec<Tvt, DGrid>,
             typename T_bd = std::array< std::array<Trl, 2>, DGrid>
             >
  bool is_migrate<T_q, T_bd> ( const T_q& q, const T_bd& bounds ) noexcept {
    auto&& tmp = apt::per_dim::make<N>
      ( []( const auto& x, const auto& bd ) noexcept {
          return x < std::get<0>(bd) || x > std::get<1>(bd);
        }, q, bounds );
    return std::apply( []( auto... args){ return (... || args); }, std::move(tmp) );
  }

  namespace impl {
    template < std::size_t DGrid, typename T, std::size_t DPtc, std::size_t I >
    inline void migrate( particle::vector<T,DPtc>& buffer,
                  const std::array< std::array<int,2>, DGrid >& neighbors,
                  const std::array< std::array<T,2>, DGrid>& bounds,
                  const mpi::Comm& comm ) {
      auto&& lcr =
        [&bd=std::get<I>(bounds)]( const auto& ptc ) noexcept {
          return ( std::get<I>(ptc.q) >= std::get<0>(bd) ) + ( std::get<I>(ptc.q) > std::get<1>(bd) );
        };
      impl::migrate_1dim( buffer, std::get<I>(neighbors), std::move(lcr), comm );
      if constexpr ( I > 0 )
                     migrate<DGrid, T, DPtc, I-1>( buffer, neighbors, bounds, comm );
    }
  }

  template < std::size_t DGrid, typename T, std::size_t DPtc,
             typename PtcVector = particle::vector<T,DPtc>,
             typename T_neigh = std::array< std::array<int,2>, DGrid >,
             typename T_bd = std::array< std::array<T,2>, DGrid>,
             typename Comm = mpi::Comm
             >
  inline void migrate<PtcVector, T_neigh, T_bd, Comm>( PtcVector& buffer, const T_neigh& neighbors, const T_bounds& bounds, const Comm& comm ) {
    return impl::migrate<DGrid, T, DPtc, DGrid - 1>( buffer, neighbors, bounds, comm );
  }


}
