#pragma once

#include "apt/array.hpp"
#include "mpipp/mpi++.hpp"

namespace dye {

// the following are ensemble specs, which will be stored on primary and be
// passed on to all replicas
template <int DGrid>
struct Ensemble {
  mpi::Comm intra;
  apt::array<apt::pair<std::optional<mpi::InterComm>>, DGrid> inter;

  int chief = 0;  // the ensemble rank of the primary process in this ensemble
  int chief_cart_rank = 0;

  // need the following so replicas can also know
  apt::array<int, DGrid> cart_coords;
  apt::array<mpi::Topo, DGrid> cart_topos;

  int label() const noexcept;
  inline auto size() const noexcept { return intra.size(); }
  inline bool is_chief() const noexcept { return intra.rank() == chief; }
  apt::pair<bool> is_at_boundary(int ith_dim) const noexcept;
  apt::array<apt::pair<bool>, DGrid> is_at_boundary() const noexcept;

  // a convinience function NOTE the reduction is in-place, which is what most
  // of use cases is
  template <typename T>
  inline void reduce_to_chief(mpi::by op, T* buffer, int count) const {
    intra.reduce<true>(op, chief, buffer, count);
  }
};

}  // namespace dye

namespace dye {

template <int DGrid>
std::optional<Ensemble<DGrid>> create_ensemble(
    const std::optional<mpi::CartComm>& cart,
    const std::optional<mpi::Comm>& intra);

template <int DGrid>
std::optional<Ensemble<DGrid>> create_ensemble(
    const std::optional<mpi::CartComm>&
        cart);  // create ensemble that only consists of chief itself

}  // namespace dye
