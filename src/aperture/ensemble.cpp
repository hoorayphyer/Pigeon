#include "./ensemble.hpp"
#include <stdexcept>

namespace aperture {
  std::optional<mpi::CartComm> create_primary_comm( std::vector<int> dims, std::vector<bool> periodic ) {
    std::optional<mpi::CartComm> result;

    const int ndims = std::min( dims.size(), periodic.size() );

    int nprmy = 1;
    for ( auto i : dims )
      nprmy *= i;

    if ( nprmy > mpi::world.size() )
      throw std::invalid_argument("Size of the Cartesian topology exceeds the size of world!");

    // simply use first few members in mpi::world as primaries
    std::vector<int> r_prmy(nprmy);
    for ( int i = 0; i < r_prmy.size(); ++i )
      r_prmy[i] = i;

    auto comm_tmp = mpi::world.split(r_prmy);

    if ( comm_tmp ) {
      result.emplace( *comm_tmp, dims, periodic );
    }

    return result;
  }


  namespace impl {
    // NOTE InterComm constructor is blocking
    template < int DGrid >
    apt::array< apt::pair<std::optional<mpi::InterComm>>, DGrid >
    link_neighbors( const std::optional<mpi::CartComm>& cart, const mpi::Comm& intra, const apt::array< apt::pair<std::optional<int>>, DGrid >& neigh_cart_ranks ) {
      // TODO
      // const auto& crk_chief = chief_cart_rank;

      // // TODO check for deadlock
      // // loop over each primary, for each of which loop over DGrid, for each of which loop over left and right.
      // // Check first whether neighbor is null, then whether there is already an intercomm.
      // for ( int i_cart = 0; i_cart < cart->size(); ++i_cart ) { // TODO cart->size()???
      //   for ( int i_dim = 0; i_dim < DGrid; ++i_dim ) {
      //     const auto& neighs = neigh_cart_ranks[i_dim];// this is in cartesian ranks
      //     auto& itc = inter[i_dim];
      //     for ( int lr = 0; lr < 2; ++lr ) {
      //       if ( crk_chief == i_cart ) {
      //         if ( !neighs[lr] || itc[lr] ) continue;
      //         itc[lr].emplace( intra, chief, cart, *neighs[lr], 147 );
      //       } else {
      //         if ( !neighs[1-lr] || *neighs[1-lr] != i_cart || itc[1-lr] ) continue;
      //         itc[1-lr].emplace( intra, chief, cart, *neighs[1-lr], 147 );
      //       }
      //     }
      //   }
      // }

      // return result;
    }
  }

  template < int DGrid >
  std::optional<Ensemble<DGrid>> create_ensemble( const std::optional<mpi::CartComm>& cart_comm, const std::optional<mpi::Comm>& intra_opt ) {
    std::optional<Ensemble<DGrid>> ens_opt;
    if ( !intra_opt ) return ens_opt;
    const auto& intra = *intra_opt;

    ens_opt.emplace();
    auto& ensemble = *ens_opt;
    ensemble.intra = intra;

    // find out who is chief
    int chief = 0;
    ensemble.chief = chief;

    constexpr int neigh_null = -147;
    apt::array< int, 3 * 1 + 2 * DGrid + 1 * (2*DGrid) > buf; // TODOL need reflection to implement the non-hard-coded version
    if ( cart_comm ) {
      const auto& cart = *cart_comm;
      ensemble.label = cart.linear_coord();
      ensemble.chief_cart_rank = cart.rank();

      auto[ coords, dims ] = cart.coords_dims();
      for ( int i = 0; i < DGrid; ++i ) {
        ensemble.cart_coords[i] = coords[i];
        ensemble.cart_dims[i] = dims[i];
        ensemble.neigh_cart_ranks[i] = cart.shift( i, 1 );
      }

      buf[0] = ensemble.label;
      buf[1] = ensemble.chief_cart_rank;
      for ( int i = 0; i < DGrid; ++i ) {
        buf[2+i] = ensemble.cart_coords[i];
        buf[2+DGrid+i] = ensemble.cart_dims[i];
        for ( int lr = 0; lr < 2; ++lr )
          buf[2+2*DGrid+2*i+lr] = ensemble.neigh_cart_ranks[i][lr] ? *(ensemble.neigh_cart_ranks[i][lr]) : neigh_null;
      }

      intra.broadcast( ensemble.chief, buf.data(), buf.size() );

    } else {
      intra.broadcast( ensemble.chief, buf.data(), buf.size() );
      ensemble.label = buf[0];
      ensemble.chief_cart_rank = buf[1];
      for ( int i = 0; i < DGrid; ++i ) {
        ensemble.cart_coords[i] = buf[2+i];
        ensemble.cart_dims[i] = buf[2+DGrid+i];
        for ( int lr = 0; lr < 2; ++lr ) {
          const auto& n = buf[2+2*DGrid+2*i+lr];
          if ( n != neigh_null ) ensemble.neigh_cart_ranks[i][lr].emplace(n);
        }
      }
    }

    ensemble.inter = impl::link_neighbors( cart_comm, intra, ensemble.neigh_cart_ranks );

    return ens_opt;
  }

  template < int DGrid >
  apt::array< apt::pair<bool>, DGrid >
  Ensemble<DGrid>::is_at_boundary() const {
    apt::array< apt::pair<bool>, DGrid > res;
    for ( int i = 0; i < DGrid; ++i ) {
      for ( int lr = 0; lr < 2; ++lr )
        res[i][lr] = !neigh_cart_ranks[i][lr];
    }
  }

}
