#include "parallel/locale.hpp"
#include <stdexcept>

namespace parallel {
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

  template < int DGrid >
  std::optional<Locale<DGrid>> create_locale( const std::optional<mpi::CartComm>& cart_comm, const std::optional<mpi::Comm>& ensemble_opt ) {
    std::optional<Locale<DGrid>> locale_opt;
    if ( !ensemble_opt ) return locale_opt;
    const auto& ensemble = *ensemble_opt;

    locale_opt.emplace();
    auto& locale = *locale_opt;
    constexpr int neigh_null = -147;
    std::array<int, 2 + DGrid + 2*DGrid> buf;
    if ( cart_comm ) {
      const auto& cart = *cart_comm;
      locale.label = cart.linear_coord();
      locale.chief_cart_rank = cart.rank();
      auto coords = cart.coords();
      for ( int i = 0; i < DGrid; ++i ) {
        locale.cart_coords[i] = coords[i];
        locale.neighbors[i] = cart.shift( i, 1 );
      }

      buf[0] = locale.label;
      buf[1] = locale.chief_cart_rank;
      for ( int i = 0; i < DGrid; ++i ) {
        buf[2+i] = locale.cart_coords[i];
        for ( int lr = 0; lr < 2; ++lr )
          buf[2+DGrid+2*i+lr] = locale.neighbors[i][lr] ? *locale.neighbors[i][lr] : neigh_null;
      }

      ensemble.broadcast( buf, locale.chief );

    } else {
      ensemble.broadcast( buf, locale.chief );
      locale.label = buf[0];
      locale.chief_cart_rank = buf[1];
      for ( int i = 0; i < DGrid; ++i ) {
        locale.cart_coords[i] = buf[2+i];
        for ( int lr = 0; lr < 2; ++lr ) {
          const auto& n = buf[2+DGrid+2*i+lr];
          if ( n != neigh_null ) locale.neighbors[i][lr].emplace(n);
        }
      }
    }

    for ( int i = 0; i < DGrid; ++i ) {
      for ( int lr = 0; lr < 2; ++lr )
        locale.is_at_boundary[i][lr] = !locale.neighbors[i][lr];
    }

    return locale_opt;
  }

  // NOTE InterComm constructor is blocking
  template < int DGrid >
  std::array< std::array< std::optional<mpi::InterComm>, 2>, DGrid >
  link_neighbors( const std::optional<mpi::CartComm>& cart,
                  const std::optional<mpi::Comm>& ensemble_opt,
                  const std::optional<Locale<DGrid>>& locale_opt ) {
    std::array< std::array< std::optional<mpi::InterComm>, 2>, DGrid > result;
    if (! ensemble_opt ) return result;
    const auto& ensemble = *ensemble_opt;
    auto& locale = *locale_opt;

    const auto& crk_chief = locale.chief_cart_rank;

    // loop over each primary, for each of which loop over DGrid, for each of which loop over left and right.
    // Check first whether neighbor is null, then whether there is already an intercomm.
    for ( int i_cart = 0; i_cart < cart->size(); ++i_cart ) {
      for ( int i_dim = 0; i_dim < DGrid; ++i_dim ) {
        const auto& neighs = locale.neighbors[i_dim];// this is in cartesian ranks
        auto& intercomms = result[i_dim];
        for ( int lr = 0; lr < 2; ++lr ) {
          if ( crk_chief == i_cart ) {
            if ( !neighs[lr] || intercomms[lr] ) continue;
            intercomms[lr].emplace( ensemble, locale.chief, cart, *neighs[lr], 147 );
          } else {
            if ( !neighs[1-lr] || *neighs[1-lr] != i_cart || intercomms[1-lr] ) continue;
            intercomms[1-lr].emplace( ensemble, locale.chief, cart, *neighs[1-lr], 147 );
          }
        }
      }
    }

    return result;
  }

}
