#include <memory>

#include "apt/ternary.hpp"
#include "apt/type_traits.hpp"
#include "mpipp/mpi++.hpp"
#include "particle/migration.hpp"

namespace particle ::impl {
// NOTE Assume there is no empty particles.
template <typename Buffer, typename F_LCR>
auto lcr_sort(Buffer& buffer, const F_LCR& lcr) noexcept {
  // RATIONALE sort buffer into [CCC...LLL...RRR...)
  // Let ec points to end of CCC.., or equivalently beginning of LLL...;
  // Let el points to end of LLL...;
  // Let rer points to the reverse end of ...RRR, i.e. all Rs will be (rer,
  // buffer.size() ) the upshot: ccc.. will be [0, ec), lll... will be [ec, el),
  // rrr... will be [el, buffer.size()), edge case safe
  int ec = 0, el = 0, rer = buffer.size() - 1;
  // REQUIRE rer always points to a value different from R
  while (rer > -1 && lcr(buffer[rer]) == apt::R) --rer;

  while (el <= rer) {
    switch (lcr(buffer[el])) {
      case apt::L:
        ++el;
        break;
      case apt::C:
        if (ec != el) buffer[el].swap(buffer[ec]);
        ++ec;
        ++el;
        break;
      case apt::R:
        buffer[el].swap(buffer[rer--]);
        while (el <= rer && lcr(buffer[rer]) == apt::R) --rer;
        break;
    }
  }

  return apt::array<const int, 3>{ec, el, buffer.size()};
}

template <typename Ptc, typename LCR_on_Ptc>
void migrate_1dim(std::vector<Ptc>& buffer,
                  const apt::pair<std::optional<mpi::InterComm>>& intercomms,
                  const LCR_on_Ptc& lcr, unsigned int shift) {
  // sort order is center | left | right | empty. Returned are the delimiters
  // between these catogories
  const auto begs =
      lcr_sort(buffer, lcr);  // begs = { begL, begR, begE_original };

  // NOTE buffer may be relocated in response to size growing. DO NOT store its
  // pointers or references
  int begE_run = begs[2];  // running begin of empty particles in buffer
  for (int lr = 0; lr < 2; ++lr) {
    const int tag = 147 + lr;
    std::vector<mpi::Request> reqs;
    // sending
    if (intercomms[lr]) {
      const auto& send_comm = *(intercomms[lr]);
      int local_rank = send_comm.rank();
      int remote_dest = (local_rank + shift) % send_comm.remote_size();
      reqs.push_back(send_comm.Isend(remote_dest, tag, buffer.data() + begs[lr],
                                     begs[lr + 1] - begs[lr]));
    }

    // receiving
    std::unique_ptr<Ptc[]> p_tmp(nullptr);
    int tot_num_recv = 0;
    if (intercomms[1 - lr]) {
      const auto& recv_comm = *(intercomms[1 - lr]);
      int local_rank = recv_comm.rank();
      int local_size = recv_comm.size();
      int remote_size = recv_comm.remote_size();
      std::vector<int> remote_srcs;
      std::vector<int> scan_recv_counts = {0};  // exclusive scan
      int src_rank =
          (local_rank + local_size - (shift % local_size)) % local_size;
      while (src_rank < remote_size) {
        remote_srcs.push_back(src_rank);
        scan_recv_counts.push_back(scan_recv_counts.back() +
                                   recv_comm.probe(src_rank, tag, buffer[0]));
        src_rank += local_size;
      }

      tot_num_recv = scan_recv_counts.back();
      // If recved more than space allows, store them in a temporary buffer then
      // later merge with the primary buffer
      Ptc* p_recv = nullptr;
      if (begE_run + tot_num_recv <= buffer.capacity()) {
        buffer.resize(begE_run + tot_num_recv);
        p_recv = buffer.data() + begE_run;
      } else {
        p_tmp.reset(new Ptc[tot_num_recv]);
        p_recv = p_tmp.get();
      }

      for (int i = 0; i < remote_srcs.size(); ++i) {
        reqs.push_back(
            recv_comm.Irecv(remote_srcs[i], tag, p_recv + scan_recv_counts[i],
                            scan_recv_counts[i + 1] - scan_recv_counts[i]));
      }
    }

    mpi::waitall(reqs);

    // merge into buffer if needed.
    if (p_tmp) {
      buffer.resize(begE_run + tot_num_recv);
      for (int i = 0; i < tot_num_recv; ++i)
        buffer[begE_run + i] = std::move(p_tmp[i]);
      p_tmp.reset(nullptr);
    }
    begE_run += tot_num_recv;
  }

  {                           // overwrite sent particles with recved ones
    int itr = begs[0];        // forward iterator
    int ritr = begE_run - 1;  // reverse iterator
    while (itr != begs[2] && ritr != begs[2] - 1) {
      std::swap(buffer[itr++], buffer[ritr--]);
    }
  }

  // NOTE It is essential to not have empty particles within buffer.size
  // otherwise lcr_sort will fail.
  buffer.resize(begs[0] + begE_run - begs[2]);
}
}  // namespace particle::impl

namespace particle {
template <typename Real, template <typename> class PtcSpecs, int DGrid,
          int I = DGrid - 1>
void migrate(std::vector<Particle<Real, PtcSpecs>>& buffer,
             const apt::array<mpi::Topo, DGrid>& topos,
             const apt::array<apt::pair<std::optional<mpi::InterComm>>, DGrid>&
                 intercomms,
             unsigned int pairing_shift) {
  // NOTE when topo[I] is a single-cpu ring, there is no going outside of bulk.
  // Hence nothing needs to be done.
  if (!(topos[I].periodic() && topos[I].dim() == 1)) {
    auto lcr = [](const auto& ptc) noexcept {
      return (ptc.template get<migrcode, DGrid>() % apt::pow3(I + 1)) /
             apt::pow3(I);
    };

    impl::migrate_1dim(buffer, intercomms[I], lcr, pairing_shift);
  }
  if constexpr (I > 0)
    migrate<Real, PtcSpecs, DGrid, I - 1>(buffer, topos, intercomms,
                                          pairing_shift);
}
}  // namespace particle
