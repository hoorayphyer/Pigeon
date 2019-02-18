#ifndef  _LOCALE_HPP_
#define  _LOCALE_HPP_

#include "mpi++.hpp"

namespace parallel {

  // the following are ensemble specs, which will be stored on primary and be passed on to all replicas
  template < int DGrid >
  struct Locale {
    const int chief = 0; // the ensemble rank of the primary process in this ensemble
    int label = 0; // ensemble label;
    int chief_cart_rank = 0;

    // need the following so replicas can also know
    std::array< int, DGrid > cart_coords;
    std::array< std::array< std::optional<int>, 2 >, DGrid > neighbors; // cartesian rank, for each dimension, 0 means upstream, 1 means downstream
    std::array< std::array<bool, 2>, DGrid > is_at_boundary;

    // std::array< std::array<bool, 2>, DGrid > is_axis;

    // std::array< int, DGrid > anchor;
    // std::array< int, DGrid > extent;
  };

}

namespace parallel {
  std::optional<mpi::Comm> create_primary_comm ( std::vector<int> dims, std::vector<bool> periodic );
  std::optional<mpi::Comm> create_ensemble_comm ( const std::vector<int>& members );

  template < int DGrid >
  Locale<DGrid> create_locale( const std::optional<mpi::Comm>& cart, const mpi::Comm& ensemble );

  template < int DGrid >
  std::array< std::array< std::optional<mpi::InterComm>, 2>, DGrid >
  link_neighbors( const std::optional<mpi::Comm>& cart,
                  const mpi::Comm& ensemble,
                  const Locale<DGrid>& locale );

}

#endif