#ifndef  _COMMUNICATOR_HPP_
#define  _COMMUNICATOR_HPP_

#include <memory>
#include <vector>
#include <optional>
#include "mpi++.hpp"

namespace parallel {
  std::optional<mpi::Comm> create_primary_comm ( std::vector<int> dims, std::vector<bool> periodic );

  // TODO link label and root using caching in mpi
  std::optional<mpi::Comm> create_ensemble_comm ( int label, int chief, const std::vector<int>& members );

  namespace ensemble {
    int chief( const mpi::Comm& ens );
    int label( const mpi::Comm& ens );
  }
}

#endif
