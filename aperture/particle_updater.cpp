#include "particle_updater.hpp"
#include "dye/ensemble.hpp"

#include "field/mesh_shape_interplay.hpp"

#include "particle_properties.hpp"
#include "particle/forces.hpp"
#include "particle/scattering.hpp"
#include "particle/migration.hpp"

#include "kernel/coordinate.hpp"

#include "kernel/grid.hpp"
#include "kernel/shapef.hpp"

namespace particle {
  map<Properties> properties;

  template < class PtcArr >
  ScatGen<PtcArr> scat_gen;

  template < class Ptc >
  ForceGen<Ptc> force_gen;
}

// TODOL all the stuff under this {} are meant to be user-specified. Here the pulsar in LogSpherical is used
namespace particle {
  template < class Ptc >
  struct CurvatureRadiate : public scat::Channel<Ptc> {
    using T = typename Ptc::vec_type::element_type;
    using Vec = apt::Vec<T,Ptc::NDim>;

    constexpr T calc_Rc ( const Ptc& ptc, const Vec& dp, T dt, const Vec& B ) noexcept {
      // qB / (\gamma m c) * dt < 2 \pi / 10
      bool is_gyration_resolved =
        std::sqrt( apt::sqabs(B) / ( (mass_x != 0) + apt::sqabs(ptc.p()) ) ) * abs_charge_x * dt / mass_x < (2 * std::acos(-1) / 10.0);

      if ( is_gyration_resolved ) {
        // find momentum at half time step
        auto phalf = ptc.p() - dp * 0.5; // NOTE ptc.p() is already the updated p
        auto v = phalf / std::sqrt( (mass_x != 0) + apt::sqabs(phalf) );
        auto vv = apt::sqabs(v);
        Vec a = dp / dt; // a is for now force, will be converted to dv/dt
        // convert a to dv/dt
        a = ( a - v * apt::dot(v,a) ) * std::sqrt( 1.0 - vv );
        auto va = apt::dot(v,a); // get the real v dot a
        return vv / std::max( std::sqrt( apt::sqabs(a) - va * va / vv ), 1e-6 ); // in case denominator becomes zero
      } else {
        // Dipolar radius of curvature in LogSpherical
        const auto& theta = ptc.q()[1];
        auto tmp = 2.5 + 1.5 * std::cos( 2 * theta );
        return std::exp(ptc.q()[0]) * tmp * std::sqrt(tmp) / ( ( tmp + 2.0 ) * std::sin(theta) );
      }
    }

    constexpr T sample_E_ph() noexcept { return 3.5; }

    unsigned int abs_charge_x{};
    unsigned int mass_x{};

    T K_thr{};
    T gamma_off{};
    T emission_rate{};

    // return sampled energy if any
    std::optional<T> operator() ( const Ptc& ptc, const apt::Vec<T,Ptc::NDim>& dp, T dt,
                                  const apt::Vec<T,Ptc::NDim>& B, util::Rng<T>& rng ) override {
      T Rc = calc_Rc( ptc, dp, dt, B );
      T gamma = std::sqrt( (mass_x != 0) + apt::sqabs(ptc.p()) );

      if (gamma > gamma_off && gamma > K_thr *  std::cbrt(Rc) && rng.uniform() < emission_rate * dt * gamma  / Rc ) {
        return { std::min( sample_E_ph(), gamma - 1.0 ) };
      } else return {};

    }
  };

  template < class Ptc >
  struct MagneticConvert : public scat::Channel<Ptc> {
    using T = typename Ptc::vec_type::element_type;
    T B_thr{};
    T mfp{};

    std::optional<T> operator() ( const Ptc& photon, const apt::Vec<T,Ptc::NDim>& dp, T dt,
                                  const apt::Vec<T,Ptc::NDim>& B, util::Rng<T>& rng ) override {
      // prob_mag_conv = dt / mfp_mag_conv
      return ( apt::sqabs(B) > B_thr * B_thr ) && ( rng.uniform() < dt / mfp )  ? std::optional<T>(0.0) : std::nullopt;
    }
  };

  template < class Ptc >
  struct TwoPhotonCollide : public scat::Channel<Ptc> {
    using T = typename Ptc::vec_type::element_type;
    // TODO double check this implementation, it is not equivalent because the original one has some sort of gaussian in it. Use Monte Carlo
    // inline T f_x ( T x ) {
    //   // distribution of x*exp(-x^2/2), which peaks at x = 1.
    //   return std::sqrt( -2.0 * std::log(x) );
    // }
    T mfp{};

    std::optional<T> operator() ( const Ptc& photon, const apt::Vec<T,Ptc::NDim>& dp, T dt,
                                  const apt::Vec<T,Ptc::NDim>& B, util::Rng<T>& rng ) override {
      // prob_mag_conv = dt / mfp_ph_ph
      return ( rng.uniform() < dt / mfp )  ? std::optional<T>(0.0) : std::nullopt;
    }
  };


  // LogSpherical
  namespace force {
    template < typename Ptc >
    void gravity( Ptc& ptc, ts::Real<Ptc> dt, const ts::Vec<Ptc>& , const ts::Vec<Ptc>& , ts::Real<Ptc> g  ) noexcept {
      ptc.p()[0] -= g * std::exp( - 2 * ptc.q()[0] ) * dt;
    };
  }
}

namespace particle {
  template < typename Real, int DGrid, int DPtc, typename state_t, typename ShapeF,
             typename Real_dJ, knl::coordsys CS >
  ParticleUpdater< Real, DGrid, DPtc, state_t, ShapeF, Real_dJ, CS >
  ::ParticleUpdater( const knl::Grid< Real, DGrid >& localgrid, const util::Rng<Real>& rng, const std::optional<mpi::CartComm>& cart, const dye::Ensemble<DGrid>& ensemble )
    : _localgrid(localgrid),
      _dJ( knl::dims(localgrid), ShapeF() ),
      _rng(rng), _cart(cart), _ensemble(ensemble) {
    using PtcArr = array<Real, DPtc, state_t>;
    using Ptc = typename PtcArr::particle_type;

    Real mu0 = 60000.0;
    {
      Real landau0_B_thr = 0.1 * mu0;
      using force_t = force::force_t<Ptc>;
      {
        auto sp = species::electron;
        force::Force<Ptc> force;
        const auto& prop = properties.at(sp);

        force.add( force::template lorentz<Ptc>, static_cast<Real>(prop.charge_x) / prop.mass_x );
        force.add( force::gravity<Ptc>, 1.8 );
        force.add( force::template landau0<Ptc>, landau0_B_thr );

        force_gen<Ptc>.Register( sp, force );
      }
      {
        auto sp = species::positron;
        force::Force<Ptc> force;
        const auto& prop = properties.at(sp);

        force.add( force::template lorentz<Ptc>, static_cast<Real>(prop.charge_x) / prop.mass_x );
        force.add( force::gravity<Ptc>, 1.8 );
        force.add( force::template landau0<Ptc>, landau0_B_thr );

        force_gen<Ptc>.Register( sp, force );
      }
      {
        auto sp = species::ion;
        force::Force<Ptc> force;
        const auto& prop = properties.at(sp);

        force.add( force::template lorentz<Ptc>, static_cast<Real>(prop.charge_x) / prop.mass_x );
        force.add( force::gravity<Ptc>, 1.8 );
        // force.add( force::template landau0<Ptc>, landau0_B_thr );

        force_gen<Ptc>.Register( sp, force );
      }
    }


    {
      scat::RadiationFromCharges<false,PtcArr> ep_scat;

      {
        struct InZone : public scat::Eligible<Ptc> {
          bool operator() ( const Ptc& ptc ) override {
            return ptc.q()[0] < std::log(9.0);
          }
        };
        ep_scat.add(InZone{});
      }

      {
        CurvatureRadiate<Ptc> cr;
        cr.abs_charge_x = 1u;
        cr.mass_x = 1u;
        cr.K_thr = 20.0;
        cr.gamma_off = 15.0;
        cr.emission_rate = 0.25;

        ep_scat.add(cr);
      }

      scat_gen<PtcArr>.Register( species::electron, ep_scat );
      scat_gen<PtcArr>.Register( species::positron, ep_scat );
    }

    {
      scat::PhotonPairProduction<PtcArr> photon_scat;
      {
        // TODOL are there any restrictions on where photons can or cannot produce pairs?
        // photon_scat.add(InZone<Ptc>{});
      }

      {
        MagneticConvert<Ptc> mc;
        mc.B_thr = 1000.0;
        mc.mfp = 0.2;
        photon_scat.add(mc);
      }
      {
        TwoPhotonCollide<Ptc> tc;
        tc.mfp = 5.0;
        photon_scat.add(tc);
      }

      scat_gen<PtcArr>.Register( species::photon, photon_scat );
    }

  }
}

namespace particle {
  template < typename Real, int DGrid, int DPtc, typename state_t, typename ShapeF,
             typename Real_dJ, knl::coordsys CS >
  void ParticleUpdater< Real, DGrid, DPtc, state_t, ShapeF, Real_dJ, CS >
  ::update_species( species sp,
                    array<Real,3,state_t>& sp_ptcs,
                    Real dt, Real unit_e,
                    const field::Field<Real,3,DGrid>& E,
                    const field::Field<Real,3,DGrid>& B,
                    const apt::array< apt::pair<Real>, DGrid >& borders
                    ) {
    if ( sp_ptcs.size() == 0 ) return;

    const auto& prop = properties.at(sp);

    using PtcArr = array<Real, DPtc, state_t>;
    using Ptc = typename PtcArr::particle_type;

    auto update_p = force_gen<Ptc>(sp);

    auto* scat = scat_gen<PtcArr>(sp);

    auto shapef = ShapeF();
    auto charge_over_dt = prop.charge_x * unit_e / dt;

    auto migrate_dir =
      []( auto q, auto lb, auto ub ) noexcept {
        return ( q >= lb ) + ( q > ub );
      };

    auto update_q =
      [is_massive=(prop.mass_x != 0)] ( auto& ptc, Real dt ) {
        auto gamma = std::sqrt( is_massive + apt::sqabs(ptc.p()) );

        if constexpr ( CS == knl::coordsys::Cartesian ) {
            // a small optimization for Cartesian
            knl::coord<CS>::geodesic_move( ptc.q(), ptc.p(), dt / gamma );
          } else {
          knl::coord<CS>::geodesic_move( ptc.q(), (ptc.p() /= gamma), dt );
          ptc.p() *= gamma;
        }
      };

    auto abs2std =
      [&grid=_localgrid]( const auto& qabs ) {
        apt::array<Real,DPtc> q_std;
        apt::foreach<0,DGrid> // NOTE DGrid instead of DPtc
          ( [](auto& q, auto q_abs, const auto& g ) noexcept {
              q = ( q_abs - g.lower() ) / g.delta();
            }, q_std, qabs, grid );
        return q_std;
      };

    for ( auto ptc : sp_ptcs ) { // TODOL sematics, check ptc is proxy
      if( ptc.is(flag::empty) ) continue;

      {
        auto q0_std = abs2std( ptc.q() );
        auto E_itpl = field::interpolate( E, q0_std, shapef );
        auto B_itpl = field::interpolate( B, q0_std, shapef );

        apt::Vec<Real,DPtc> dp = -ptc.p();
        update_p( ptc, dt, E_itpl, B_itpl );
        if ( scat ) {
          dp += ptc.p();
          (*scat)( std::back_inserter(sp_ptcs), ptc, std::move(dp), dt, B_itpl, _rng );
        }

        // NOTE q is updated, starting from here, particles may be in the guard cells.
        update_q( ptc, dt );
        // TODO pusher handle boundary condition. Is it needed?
        if ( prop.charge_x != 0 )
          _dJ.deposit( charge_over_dt, q0_std, abs2std(ptc.q()) );
      }

      {
        char mig_dir = 0;
        for ( int i = 0; i < DGrid; ++i ) {
          mig_dir += migrate_dir( ptc.q()[i], borders[i][LFT], borders[i][RGT] ) * pow3(i);
        }

        if ( mig_dir != (pow3(DGrid) - 1) / 2 ) {
          _migrators.emplace_back(std::move(ptc));
          _migrators.back().extra() = mig_dir;
        }
      }

    }
  }
}

namespace particle {
  template < typename Real, int DGrid, int DPtc, typename state_t,
             typename ShapeF,
             typename Real_dJ,
             knl::coordsys CS
             >
  void ParticleUpdater< Real, DGrid, DPtc, state_t, ShapeF, Real_dJ, CS >
  ::operator() ( field::Field<Real,3,DGrid>& J,
                 map<array<Real,3,state_t>>& particles,
                 const field::Field<Real,3,DGrid>& E,
                 const field::Field<Real,3,DGrid>& B,
                 const apt::array< apt::pair<Real>, DGrid >& borders,
                 Real dt, Real unit_e, int timestep ) {

    _dJ.reset();

    // TODO NOTE injection and annihilation are more like free sources and sinks of particles users control, in other words they fit the notion of particle boundary conditions. Do them outside of Updater. In contrast, pair creation is a physical process we model, so it is done inside the updater.
    // inject_particles(); // TODO only coordinate space current is needed in implementing current regulated injection // TODO compatible with ensemble??
    // TODOL annihilation will affect deposition // NOTE one can deposit in the end
    // annihilate_mark_pairs( );

    for ( auto&[ sp, ptcs ] : particles ) {
      const auto old_size = ptcs.size();
      update_species( sp, ptcs, dt, unit_e, E, B, borders );

      // Put particles where they belong after scattering
      for ( auto i = old_size; i < ptcs.size(); ++i ) {
        auto this_sp = ptcs[i].template get<species>();
        if ( this_sp != sp ) {
          particles[this_sp].push_back(std::move(ptcs[i]));
        }
      }
    }

    migrate( _migrators, _ensemble.inter, timestep );
    for ( auto&& ptc : _migrators ) {
      particles[ptc.template get<species>()].push_back( std::move(ptc) );
    }
    _migrators.resize(0);

    _dJ.reduce( _ensemble.chief, _ensemble.intra );

    if ( _cart ) {
      const auto& Jmesh = _dJ.integrate( *_cart );

      { // get J from Jmesh NOTE one needs to rescale Jmesh back to real grid delta
        const auto& mesh = J.mesh();
        const auto& grid = _localgrid;
        using coord_t = knl::coord<CS>;

        for ( int iJ = 0; iJ < 3; ++iJ ) {
          auto Jcomp = J[iJ]; // TODOL semantics;
          auto Jmesh_comp = Jmesh[iJ]; // TODOL semantics;
          const auto& offset = Jcomp.offset();

          Real delta = iJ < DGrid ? grid[iJ].delta() : 1.0;

          const decltype( coord_t::template h<0,Real> ) * scale = nullptr;
          switch( iJ < DGrid ? iJ : -1 ) {
          case 0 : scale = coord_t::template hh<0,Real>; break;
          case 1 : scale = coord_t::template hh<1,Real>; break;
          case 2 : scale = coord_t::template hh<2,Real>; break;
          case -1 : scale = coord_t::template hhh<Real>; break;
          default: scale = nullptr;
          }

          apt::array<Real,3> qs{};
          for( const auto& I : apt::Block( mesh.bulk_extent() ) ) {
            // TODOL also use iterator for qs generation
            for ( int i = 0; i < DGrid; ++i )
               qs[i] = grid[i].absc(I[i], offset[i]);

            Jcomp(I) = Jmesh_comp(I) * delta / scale( qs[0], qs[1], qs[2] );
          }
        }
      }

    }
  }

}

#include "traits.hpp"
using namespace traits;
namespace particle {
  template class ParticleUpdater< real_t, DGrid, DPtc, ptc_state_t, ShapeF,
                                  real_dj_t, coordinate_system >;
}
