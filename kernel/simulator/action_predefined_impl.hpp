#include <vector>

#include "particle/array.hpp"
#include "particle/migration.hpp"
#include "particle/properties.hpp"
#include "simulator/action_predefined.hpp"

namespace particle {

template <int DGrid, typename R, template <typename> class S, typename RJ>
void Migrator<DGrid, R, S, RJ>::impl(map<array<R, S>>& particles,
                                     std::vector<Particle<R, S>>& new_ptc_buf,
                                     const map<Properties>& properties,
                                     const apt::Grid<R, DGrid>& grid,
                                     const dye::Ensemble<DGrid>& ens,
                                     int timestep) const {
  // bulk range = [lb, ub)
  constexpr auto migrtrit = [](auto q, auto lb, auto ub) noexcept {
    return (q >= lb) + (q >= ub);
  };

  for (auto sp : particles) {
    for (auto ptc : particles[sp]) {  // TODOL semantics
      if (!ptc.is(flag::exist)) continue;
      int mig_dir{};
      for (int i = 0; i < DGrid; ++i) {
        mig_dir +=
            migrtrit(ptc.q(i), grid[i].lower(), grid[i].upper()) * apt::pow3(i);
      }

      if (mig_dir != (apt::pow3(DGrid) - 1) / 2) {
        ptc.template set<migrcode, DGrid>(mig_dir);
        new_ptc_buf.emplace_back(std::move(ptc));
      }
    }
  }

  migrate(new_ptc_buf, ens.cart_topos, ens.inter, timestep);

  // NOTE adjust particle positions in the ring topology, regardless how many
  // cpus there are on that ring.
  {
    bool has_periodic = false;
    for (int i = 0; i < DGrid; ++i) {
      has_periodic = has_periodic || ens.cart_topos[i].periodic();
    }
    if (has_periodic) {
      for (auto& ptc : new_ptc_buf) {
        if (!ptc.is(flag::exist)) continue;

        for (int i = 0; i < DGrid; ++i) {
          if (!ens.cart_topos[i].periodic()) continue;
          int idx = static_cast<int>((ptc.q(i) - m_supergrid[i].lower()) /
                                         m_supergrid[i].delta() +
                                     0.5);
          if (idx >= 0)
            idx /= m_supergrid[i].dim();
          else
            idx = -((-idx) / m_supergrid[i].dim() + 1);

          ptc.q(i) -= idx * m_supergrid[i].dim() * m_supergrid[i].delta();
        }
      }
    }
  }

  for (auto&& ptc : new_ptc_buf) {
    if (!ptc.is(flag::exist)) continue;
    auto sp = ptc.template get<species>();
#if PIC_DEBUG
    // // check if the received ptc trully resides in this ensemble.
    // apt::array<int,DGrid> mig_co;
    // bool is_OK = true;
    // for ( int i = 0; i < DGrid; ++i ) {
    //   mig_co[i] = migrate_code( ptc.q(i), _borders[i][LFT],
    //   _borders[i][RGT] ); if ( mig_co[i] != 1 &&
    //   !_ens_opt->is_at_boundary(i)[(mig_co[i] != 0)] ) //NOTE need to
    //   consider boundaries
    //     is_OK = false;
    // }
    // if ( !is_OK ) {
    //   lgr::file << "ts=" << debug::timestep << ", wr=" << debug::world_rank
    //   << ", el=" << debug::ens_label << std::endl; lgr::file << "Received
    //   across-ensemble particles! q = " << ptc.q() << ", p = " << ptc.p() <<
    //   std::endl; lgr::file << "  mig_dir on new ensemble  = " << mig_co;
    //   // get old mig_co
    //   for ( int i = 0; i < DGrid; ++i ) {
    //     mig_co[i] = ( migrInt<DGrid>(ptc) % apt::pow3(i+1) ) /
    //     apt::pow3(i);
    //   }
    //   lgr::file << ", mig_dir on old ensemble = " << mig_co << std::endl;
    //   debug::throw_error("Received across-ensemble particles!");
    // }
#endif
    ptc.template reset<migrcode>();
    particles[sp].push_back(std::move(ptc));
  }
  new_ptc_buf.resize(0);
}
}  // namespace particle
