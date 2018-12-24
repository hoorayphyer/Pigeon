#include "parallel/communicator.hpp"
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
      mpi::topo::cartesianize( *result, dims, periodic );
    }

    return result;
  }

  std::optional<mpi::Comm> create_ensemble_comm ( int label, int chief, const std::vector<int>& members ) {
    std::optional<mpi::Comm> result;

    if ( std::find( members.begin(), members.end(), mpi::world.rank() ) != members.end() ) {
      result.emplace( members );
      (*result).attrs["label"] = label;
      (*result).attrs["chief"] = chief;
    }
    return result;
  }

  namespace ensemble {
    int chief(mpi::Comm comm) {
      return std::any_cast<int>(comm.attrs.at("chief"));
    }

    int label(mpi::Comm comm) {
      return std::any_cast<int>(comm.attrs.at("label"));
    }
  }
}
