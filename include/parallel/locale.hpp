#ifndef  _LOCALE_HPP_
#define  _LOCALE_HPP_

#include "mpi++.hpp"
#include "apt/array.hpp"

namespace parallel {

  // the following are ensemble specs, which will be stored on primary and be passed on to all replicas
  template < int DGrid >
  struct Locale {
    const int chief = 0; // the ensemble rank of the primary process in this ensemble
    int label = 0; // ensemble label;
    int chief_cart_rank = 0;

    // need the following so replicas can also know
    apt::array< int, DGrid > cart_coords;
    apt::array< apt::pair<std::optional<int>>, DGrid > neighbors; // cartesian rank, for each dimension, 0 means upstream, 1 means downstream
    apt::array< apt::pair<bool>, DGrid > is_at_boundary;

  };

}

namespace parallel {
  std::optional<mpi::CartComm> create_primary_comm ( std::vector<int> dims, std::vector<bool> periodic );

  template < int DGrid >
  std::optional<Locale<DGrid>> create_locale( const std::optional<mpi::CartComm>& cart, const std::optional<mpi::Comm>& ensemble );

  template < int DGrid >
  apt::array< apt::pair<std::optional<mpi::InterComm>>, DGrid >
  link_neighbors( const std::optional<mpi::CartComm>& cart_comm,
                  const std::optional<mpi::Comm>& ensemble,
                  const std::optional<Locale<DGrid>>& locale );

}

#endif
