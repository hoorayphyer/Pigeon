#include "parallel/locale.hpp"
#include <stdexcept>
#include <algorithm> // std::find

namespace parallel {
  std::optional<mpi::Comm> create_primary_comm( std::vector<int> dims, std::vector<bool> periodic ) {
    std::optional<mpi::Comm> result;

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

    if ( std::find( r_prmy.begin(), r_prmy.end(), mpi::world.rank() ) != r_prmy.end() ) {
      result.emplace( r_prmy );
      mpi::cartesianize( *result, dims, periodic );
    }

    return result;
  }

  std::optional<mpi::Comm> create_ensemble_comm ( const std::vector<int>& members ) {
    std::optional<mpi::Comm> result;

    if ( std::find( members.begin(), members.end(), mpi::world.rank() ) != members.end() )
      result.emplace( members );
    return result;
  }

  template < int DGrid >
  Locale<DGrid> create_locale( const std::optional<mpi::Comm>& cart, const mpi::Comm& ensemble ) {
    constexpr int neigh_null = -147;
    Locale<DGrid> locale;

    std::array<int, 2 + DGrid + 2*DGrid> buf;
    if ( cart ) {
      locale.label = mpi::cart::linear_coord( *cart );
      locale.chief_cart_rank = *cart.rank();
      locale.cart_coords = mpi::cart::coords( *cart );
      for ( int i = 0; i < DGrid; ++i ) {
        locale.neighbors[i] = mpi::cart::shift( *cart, i, 1 );
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

    return locale;
  }

  // NOTE InterComm constructor is blocking
  template < int DGrid >
  std::array< std::array< std::optional<mpi::InterComm>, 2>, DGrid >
  link_neighbors( const std::optional<mpi::Comm>& cart,
                  const mpi::Comm& ensemble,
                  const Locale<DGrid>& locale ) {
    std::array< std::array< std::optional<mpi::InterComm>, 2>, DGrid > result;
    const auto& crk_chief = locale.chief_cart_rank;

    // loop over each primary, for each of which loop over DGrid, for each of which loop over left and right.
    // Check first if neighbor is null, then if there is already an intercomm.
    for ( int i_cart = 0; i_cart < *cart.size(); ++i_cart ) {
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
