#ifndef  _ENSEMBLE_HPP_
#define  _ENSEMBLE_HPP_

#include "parallel/mpi++.hpp"
#include "apt/array.hpp"

namespace aperture {

  // the following are ensemble specs, which will be stored on primary and be passed on to all replicas
  template < int DGrid >
  struct Ensemble {
    mpi::Comm intra;
    apt::array< apt::pair<std::optional<mpi::InterComm>>, DGrid > inter;

    int chief = 0; // the ensemble rank of the primary process in this ensemble
    int chief_cart_rank = 0;

    // need the following so replicas can also know
    apt::array< int, DGrid > cart_coords;
    apt::array< int, DGrid > cart_dims;
    apt::array< bool, DGrid > is_periodic;

    int label() const noexcept;
    apt::pair<bool> is_at_boundary( int ith_dim ) const noexcept;
    apt::array< apt::pair<bool>, DGrid > is_at_boundary() const noexcept;
  };

}

namespace aperture {

  template < int DGrid >
  std::optional<Ensemble<DGrid>> create_ensemble( const std::optional<mpi::CartComm>& cart, const std::optional<mpi::Comm>& intra );

  template < int DGrid >
  std::optional<Ensemble<DGrid>> create_ensemble( const std::optional<mpi::CartComm>& cart ); // create ensemble only consists of chief itself

}

#endif
