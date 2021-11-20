#include "msh/mesh_shape_interplay.hpp"
#include "particle/forces.hpp"
#include "particle/scattering.hpp"
#include "particle/updater.hpp"

#if PIC_DEBUG
#include "debug/debugger.hpp"
#include "logger/logger.hpp"
#endif

#ifdef LORENTZ
#include "logger/logger.hpp"
#endif

namespace particle {
template <int DGrid, typename R, template <typename> class S, typename ShapeF,
          typename RJ>
void Updater<DGrid, R, S, ShapeF, RJ>::operator()(
    map<array<R, S>>& particles, field::Field<RJ, 3, DGrid>& J,
    std::vector<Particle<R, S>>* new_ptc_buf, const map<Properties>& properties,
    const field::Field<R, 3, DGrid>& E, const field::Field<R, 3, DGrid>& B,
    const apt::Grid<R, DGrid>& grid, R dt, int timestep, util::Rng<R>& rng) const {
  auto update_species = [&](species sp) {
    auto& sp_ptcs = particles[sp];
    if (sp_ptcs.size() == 0) return;

    const auto& prop = properties[sp];

    using Ptc = typename array<R, S>::particle_type;

    auto update_p = [forces = Force<R, S>::Get(sp)](auto& ptc, auto&&... args) {
      for (int i = 0; i < forces.forces.size(); ++i) {
        (forces.forces[i])(ptc, std::forward<decltype(args)>(args)...,
                           forces.params[i]);
      }
    };

    auto scatter = [scat = Scat<R, S>::Get(sp)](auto& ptc, auto& buf,
                                                auto&&... args) {
      bool is_eligible = true;
      for (const auto& elig : scat.eligs) {
        if (!(*elig)(ptc)) {
          is_eligible = false;
          break;
        }
      }
      if (is_eligible) {
        for (const auto& chnl : scat.channels) {
          if (auto param =
                  (*chnl)(ptc, std::forward<decltype(args)>(args)...)) {
            scat.impl(std::back_inserter(buf), ptc, *param);
            break;
          }
        }
      }
    };

    constexpr auto shapef = ShapeF();
    auto charge_over_dt = prop.charge_x / dt;
    bool is_deposit_current =
        (_ignore_current.find(sp) == _ignore_current.end() and
         std::abs(prop.charge_x) > 0.01);

    for (auto ptc : sp_ptcs) {  // TODOL sematics, check ptc is proxy
      if (!ptc.is(flag::exist)) continue;
#if PIC_DEBUG
      apt::Vec<R, S<R>::Dim> q_prev = ptc.q();
      apt::Vec<R, S<R>::Dim> p_prev = ptc.p();

      auto show_ptc = [&](const auto& ptc) {
        lgr::file << "ptc.q() = " << ptc.q() << ", ptc.p() = " << ptc.p()
                  << std::endl;
        lgr::file << "frac = " << ptc.frac() << std::endl;
        lgr::file << "sp = " << static_cast<int>(ptc.template get<species>())
                  << ", is_sec = " << ptc.is(flag::secondary) << std::endl;

        lgr::file << "ptc q_prev = " << q_prev << ", ptc p_prev = " << p_prev
                  << std::endl;
      };

      if (debug::has_nan(ptc)) {
        lgr::file << "ts=" << debug::timestep << ", wr=" << debug::world_rank
                  << ", el=" << debug::ens_label << std::endl;
        lgr::file << "NANBEGINNING, code=" << debug::has_nan(ptc) << std::endl;
        show_ptc(ptc);
        debug::throw_error("NANBEGINNING");
      }
#endif

      auto q0_std = msh::to_standard(grid, ptc.q());

      if (not ptc.is(flag::ignore_force)) {
        auto E_itpl = msh::interpolate(E, q0_std, shapef);
        auto B_itpl = msh::interpolate(B, q0_std, shapef);

        apt::Vec<R, S<R>::Dim> dp = -ptc.p();
        update_p(ptc, dt, E_itpl, B_itpl);
#if PIC_DEBUG
        if (debug::has_nan(ptc)) {
          lgr::file << "ts=" << debug::timestep << ", wr=" << debug::world_rank
                    << ", el=" << debug::ens_label << std::endl;
          lgr::file << "NANUPDATEP, code=" << debug::has_nan(ptc) << std::endl;
          show_ptc(ptc);
          lgr::file << "E_itpl = " << E_itpl << ", B_itpl = " << B_itpl
                    << std::endl;
#ifdef LORENTZ
          lgr::file << force::ostr.str() << std::endl;
#endif
          debug::throw_error("NANUPDATEP");
        } else {
          q_prev = ptc.q();
          p_prev = ptc.p();
        }
#endif

        dp += ptc.p();  // ptc.p() is the updated one but still based on the
                        // same location
        scatter(ptc, *new_ptc_buf, prop, dp, dt, B_itpl,
                rng);  // FIXME new_ptc_buf may be null?
#if PIC_DEBUG
        if (debug::has_nan(ptc)) {
          lgr::file << "ts=" << debug::timestep << ", wr=" << debug::world_rank
                    << ", el=" << debug::ens_label << std::endl;
          lgr::file << "NANSCAT, code=" << debug::has_nan(ptc) << std::endl;
          show_ptc(ptc);
          lgr::file << "E_itpl = " << E_itpl << ", B_itpl = " << B_itpl
                    << std::endl;
          debug::throw_error("NANSCAT");
        } else {
          q_prev = ptc.q();
          p_prev = ptc.p();
        }
#endif
      }

      // NOTE q is updated, starting from here, particles may be in the guard
      // cells.
      auto dq = _update_q(ptc.q(), ptc.p(), dt, (prop.mass_x > 0.01));
#if PIC_DEBUG
      if (debug::has_nan(ptc)) {
        lgr::file << "ts=" << debug::timestep << ", wr=" << debug::world_rank
                  << ", el=" << debug::ens_label << std::endl;
        lgr::file << "NANUPDATEQ, code=" << debug::has_nan(ptc) << std::endl;
        show_ptc(ptc);
        debug::throw_error("NANUPDATEQ");
      } else {
        q_prev = ptc.q();
        p_prev = ptc.p();
      }
#endif
      // FIXME pusher handle boundary condition. Is it needed?
      if (is_deposit_current) {
        // change dq to q1 in the coordinate space. NOTE it is not necessarily
        // the same as ptc.q()
        for (int i = 0; i < DGrid; ++i) {
          dq[i] /= grid[i].delta();
          dq[i] += q0_std[i];
        }
        for (int i = DGrid; i < S<R>::Dim; ++i) dq[i] += q0_std[i];
        msh::deposit(J, ptc.frac() * charge_over_dt, shapef, q0_std,
                     std::move(dq));
      }
    }
  };

  for (auto sp : particles) {
    update_species(sp);
  }

  {  // NOTE rescale Jmesh back to real grid delta
    auto dV = apt::dV(grid);

    for (int i = 0; i < DGrid; ++i) {
      R tmp = grid[i].delta() / dV;
      for (auto& elm : J[i].data()) elm *= tmp;
    }

    for (int i = DGrid; i < 3; ++i) {
      for (auto& elm : J[i].data()) elm /= dV;
    }
  }
}

}  // namespace particle
