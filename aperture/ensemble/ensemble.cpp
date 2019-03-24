#include "./ensemble.hpp"
#include <stdexcept>

namespace impl {
  apt::pair< std::optional<mpi::InterComm> >
  link_neighbors( const mpi::Comm& intra,
                  const std::optional<mpi::CartComm>& cart,
                  int ith_dim, int chief,
                  int cart_coord, int cart_dim, bool is_periodic ) {
    // NOTE InterComm constructor is blocking.
    // Edge case: cart_dim == 1; TODO what if is_periodic? This will need a whole new type of communicator to do the work
    if ( cart_dim == 1 && !is_periodic ) return {};

    apt::pair< std::optional<mpi::InterComm> > inter_comms;

    auto remote_leader =
      [&cart, ith_dim]( bool to_right ) {
        if( !cart ) return -1; // not significant
        return *(cart->shift( ith_dim, 1 ) [to_right]);
      };

    // odd-even alternate
    if ( cart_coord != cart_dim - 1 )
      inter_comms[RGT].emplace( intra, chief, cart, remote_leader(RGT), 147 );

    if ( cart_coord != 0 )
      inter_comms[LFT].emplace( intra, chief, cart, remote_leader(LFT), 741 );

    if ( is_periodic ) {
      if ( cart_coord == 0 )
        inter_comms[LFT].emplace( intra, chief, cart, remote_leader(LFT), 741 );
      else if ( cart_coord == cart_dim - 1 )
        inter_comms[RGT].emplace( intra, chief, cart, remote_leader(RGT), 147 );
    }

    return inter_comms;
  }
}

// create_ensemble
namespace aperture {
  template < int DGrid >
  std::optional<Ensemble<DGrid>> create_ensemble( const std::optional<mpi::CartComm>& cart_comm, const std::optional<mpi::Comm>& intra_opt ) {
    std::optional<Ensemble<DGrid>> ens_opt;
    if ( !intra_opt ) return ens_opt;
    const auto& intra = *intra_opt;

    ens_opt.emplace();
    auto& ens = *ens_opt;
    ens.intra = intra;

    // find out who is the chief
    {
      auto comm = intra.split( {!cart_comm}, intra.rank() );
      if ( cart_comm ) ens.chief = intra.rank();
      else {
        // deduce by elimination
        auto no_cart = comm->group().translate_all( intra.group() );
        ens.chief = intra.size() * ( intra.size() - 1 ) / 2;
        for ( auto x : no_cart ) ens.chief -= x;
      }
    }


    apt::array< int, 1 + 3 * DGrid > buf; // TODOL need reflection to implement the non-hard-coded version
    if ( cart_comm ) {
      const auto& cart = *cart_comm;
      ens.chief_cart_rank = cart.rank();

      {
        auto[c, d, p] = cart.coords_dims_periodic();
        for ( int i = 0; i < DGrid; ++i ) {
          ens.cart_coords[i] = c[i];
          ens.cart_dims[i] = d[i];
          ens.is_periodic[i] = p[i];
        }
      }

      buf[0] = ens.chief_cart_rank;
      for ( int i = 0; i < DGrid; ++i ) {
        buf[1+i] = ens.cart_coords[i];
        buf[1+DGrid+i] = ens.cart_dims[i];
        buf[1+ 2*DGrid+i] = ens.is_periodic[i];
      }
    }

    intra.broadcast( ens.chief, buf.data(), buf.size() );

    if ( !cart_comm ) {
      ens.chief_cart_rank = buf[0];
      for ( int i = 0; i < DGrid; ++i ) {
        ens.cart_coords[i] = buf[1+i];
        ens.cart_dims[i] = buf[1+DGrid+i];
        ens.is_periodic[i] = buf[1+2*DGrid+i];
      }
    }

    for ( int i = 0; i < DGrid; ++i )
      ens.inter[i] = impl::link_neighbors( intra, cart_comm, i, ens.chief, ens.cart_coords[i], ens.cart_dims[i], ens.is_periodic[i] );

    return ens_opt;
  }

  template < int DGrid >
  std::optional<Ensemble<DGrid>> create_ensemble( const std::optional<mpi::CartComm>& cart ) {
    if ( !cart ) return {};
    std::optional<mpi::Comm> intra_opt(mpi::self);
    return create_ensemble<DGrid>( cart, intra_opt );
  }
}

// accessors
namespace aperture {
  template < int DGrid >
  int Ensemble<DGrid>::label() const noexcept {
    int i = DGrid - 1;
    int result = cart_coords[i];
    for ( --i; i > -1; --i )
      result = cart_coords[i] + result * cart_dims[i];

    return result;
  }


  template < int DGrid >
  apt::pair<bool> Ensemble<DGrid>::is_at_boundary( int ith_dim ) const noexcept {
    return { !(inter[ith_dim][LFT]), !(inter[ith_dim][RGT]) };
  }

  template < int DGrid >
  apt::array< apt::pair<bool>, DGrid >
  Ensemble<DGrid>::is_at_boundary() const noexcept {
    apt::array< apt::pair<bool>, DGrid > res;
    for ( int i = 0; i < DGrid; ++i )
      res[i] = is_at_boundary(i);
  }

}

#include "traits.hpp"
namespace aperture {
  using namespace traits;
  template
  std::optional<Ensemble<DGrid>> create_ensemble<DGrid>( const std::optional<mpi::CartComm>& cart, const std::optional<mpi::Comm>& intra );
}
