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
    int label = 0; // ensemble label;
    int chief_cart_rank = 0;

    // need the following so replicas can also know
    apt::array< int, DGrid > cart_coords;
    apt::array< int, DGrid > cart_dims;

    apt::array< apt::pair<std::optional<int>>, DGrid > neigh_cart_ranks; // TODOL currently used in old_field_solver, and link_neighbors

    apt::array< apt::pair<bool>, DGrid > is_at_boundary() const;

  };

}

namespace aperture {
  std::optional<mpi::CartComm> create_primary_comm ( std::vector<int> dims, std::vector<bool> periodic );

  template < int DGrid >
  std::optional<Ensemble<DGrid>> create_ensemble( const std::optional<mpi::CartComm>& cart, const std::optional<mpi::Comm>& intra );

}

#endif
