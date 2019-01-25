#include "particle_module/migration.hpp"

namespace particle :: impl {
  template < std::size_t I, typename T, std::size_t DPtc >
  auto migsort( particle::vector<T,DPtc>& buffer, const std::array<T, 2>& bd ) noexcept {
    
  }



  template < std::size_t I, typename T, std::size_t DPtc >
  void migrate_1dim ( particle::vector<T,DPtc>& buffer,
                     const std::array<int, 2>& neigh,
                     const std::array<T, 2>& bound,
                     auto comm ) {
    const auto& nL = std::get<0>(neigh);
    const auto& nR = std::get<1>(neigh);
    // sort order is center | left | right | empty. Returned are the delimiters between these catogories
    auto[begL, begR, begE] = sort<I>( buffer, compare(bound) ); // TODO finish the sort
    auto* ptr = buffer.data();
    int num_recv = 0;

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
  template < typename Tvt, std::size_t DGrid, typename Trl = vec::remove_cvref_t<Tvt> >
  bool is_migrate( const Vec<Tvt, DGrid>& q, const std::array< std::array<Trl, 2>, DGrid>& bounds ) noexcept {
    auto&& tmp = vec::per_dim::make<N>
      ( []( const auto& x, const auto& bd ) noexcept {
          return x < std::get<0>(bd) || x > std::get<1>(bd);
        }, q, bounds );
    return std::apply( []( auto... args){ return (... || args); }, std::move(tmp) );
  }

  template < std::size_t DGrid, typename T, std::size_t DPtc, std::size_t I >
  void migrate( particle::vector<T,DPtc>& buffer,
                const std::array< std::array<int,2>, DGrid >& neighbors,
                const std::array< std::array<T,2>, DGrid>& bounds,
                auto comm ) {
    impl::migrate_1dim<I, T, DPtc>( buffer, std::get<I>(neighbors), std::get<I>(bounds), comm );
    if constexpr ( I > 0 )
      migrate<DGrid, T, DPtc, I-1>( buffer, neighbors, bounds, comm );
  }

}
