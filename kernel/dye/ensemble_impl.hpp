#include <stdexcept>

#include "dye/ensemble.hpp"

namespace dye::impl {
apt::pair<std::optional<mpi::InterComm>> link_neighbors(
    const mpi::Comm& intra, const std::optional<mpi::CartComm>& cart,
    int ith_dim, int chief, int cart_coord, int cart_dim, bool is_periodic) {
  // NOTE InterComm constructor is blocking.
  apt::pair<std::optional<mpi::InterComm>> inter_comms;

  // Edge case: cart_dim == 1. If nonperiodic, nothing needs to be created. If
  // periodic, since MPI intercommunicator doesn't allow overlapping of local
  // and remote groups, just return no neighbors for now and make applications
  // check for this case in practise. // TODO periodicity
  if (cart_dim == 1) return inter_comms;

  auto remote_leader = [&cart, ith_dim](bool to_right) {
    if (!cart) return -1;  // not significant
    return *(cart->shift(ith_dim, 1)[to_right]);
  };

  // odd-even alternate.
  // Round One. Even to right except for the last, Odd to left
  if (cart_coord % 2 == 0) {
    if (cart_coord != cart_dim - 1)
      inter_comms[RGT].emplace(intra, chief, cart, remote_leader(RGT), 147);
  } else {
    inter_comms[LFT].emplace(intra, chief, cart, remote_leader(LFT), 147);
  }

  // Round Two. Even to left except for 0, Odd to right except for the last

  if (cart_coord % 2 == 0) {
    if (cart_coord != 0)
      inter_comms[LFT].emplace(intra, chief, cart, remote_leader(LFT), 741);
  } else {
    if (cart_coord != cart_dim - 1)
      inter_comms[RGT].emplace(intra, chief, cart, remote_leader(RGT), 741);
  }

  // Connect head and tail if periodic
  if (is_periodic) {
    if (cart_coord == 0)
      inter_comms[LFT].emplace(intra, chief, cart, remote_leader(LFT), 137);
    else if (cart_coord == cart_dim - 1)
      inter_comms[RGT].emplace(intra, chief, cart, remote_leader(RGT), 137);
  }

  return inter_comms;
}
}  // namespace dye::impl

// create_ensemble
namespace dye {
template <int DGrid>
std::optional<Ensemble<DGrid>> create_ensemble(
    const std::optional<mpi::CartComm>& cart_comm,
    const std::optional<mpi::Comm>& intra_opt) {
  std::optional<Ensemble<DGrid>> ens_opt;
  if (!intra_opt) return ens_opt;
  const auto& intra = *intra_opt;

  ens_opt.emplace();
  auto& ens = *ens_opt;
  ens.intra = intra;

  // find out who is the chief
  {
    auto comm = intra.split({!cart_comm}, intra.rank());
    if (cart_comm)
      ens.chief = intra.rank();
    else {
      // deduce by elimination
      auto no_cart = comm->group().translate_all(intra.group());
      ens.chief = intra.size() * (intra.size() - 1) / 2;
      for (auto x : no_cart) ens.chief -= x;
    }
  }

  apt::array<int, 1 + 3 * DGrid>
      buf;  // TODOL need reflection to implement the non-hard-coded version
  if (cart_comm) {
    const auto& cart = *cart_comm;
    ens.chief_cart_rank = cart.rank();

    {
      auto [c, t] = cart.coords_topos();
      for (int i = 0; i < DGrid; ++i) {
        ens.cart_coords[i] = c[i];
        ens.cart_topos[i] = t[i];
      }
    }

    buf[0] = ens.chief_cart_rank;
    for (int i = 0; i < DGrid; ++i) {
      buf[1 + i] = ens.cart_coords[i];
      buf[1 + DGrid + i] = ens.cart_topos[i].dim();
      buf[1 + 2 * DGrid + i] = ens.cart_topos[i].periodic();
    }
  }

  intra.broadcast(ens.chief, buf.begin(), buf.size());

  if (!cart_comm) {
    ens.chief_cart_rank = buf[0];
    for (int i = 0; i < DGrid; ++i) {
      ens.cart_coords[i] = buf[1 + i];
      ens.cart_topos[i] = mpi::Topo(buf[1 + DGrid + i],
                                    static_cast<bool>(buf[1 + 2 * DGrid + i]));
    }
  }

  for (int i = 0; i < DGrid; ++i)
    ens.inter[i] = impl::link_neighbors(
        intra, cart_comm, i, ens.chief, ens.cart_coords[i],
        ens.cart_topos[i].dim(), ens.cart_topos[i].periodic());

  return ens_opt;
}

template <int DGrid>
std::optional<Ensemble<DGrid>> create_ensemble(
    const std::optional<mpi::CartComm>& cart) {
  if (!cart) return {};
  std::optional<mpi::Comm> intra_opt(mpi::self);
  return create_ensemble<DGrid>(cart, intra_opt);
}
}  // namespace dye

// accessors
namespace dye {
template <int DGrid>
int Ensemble<DGrid>::label() const noexcept {
  return chief_cart_rank;
}

template <int DGrid>
apt::pair<bool> Ensemble<DGrid>::is_at_boundary(int ith_dim) const noexcept {
  if (cart_topos[ith_dim].periodic())
    return {};
  else
    return {!(inter[ith_dim][LFT]), !(inter[ith_dim][RGT])};
}

template <int DGrid>
apt::array<apt::pair<bool>, DGrid> Ensemble<DGrid>::is_at_boundary()
    const noexcept {
  apt::array<apt::pair<bool>, DGrid> res;
  for (int i = 0; i < DGrid; ++i) res[i] = is_at_boundary(i);
  return res;
}

}  // namespace dye
