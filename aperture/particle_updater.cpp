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
}

namespace particle {
  template < class Ptc >
  struct InZone : public scat::Eligible<Ptc> {
    bool operator() ( const Ptc& ptc ) override {
      return ptc.q()[0] < std::log(9.0);
      // return _pane.in_productive_zone( ptc.q );
    }
  };

  template < class Ptc >
  struct CurvatureRadiate : public scat::Channel<Ptc> {
    using T = typename Ptc::vec_type::element_type;

    bool is_massive;
    T gamma_off = 15.0;
    T K_thr = 20;
    // prob_fiducial = K_curv_em_rate * dt
    T prob_fiducial = 0.1;

    // float ParticlePusher::CalculateRc( Scalar dt, const Vec3<MOM_TYPE> &p, const Vec3<MOM_TYPE> &dp) const {
    //   // find momentum at half time step
    //   Vec3<MOM_TYPE> phalf( p + dp * 0.5 );
    //   Vec3<MOM_TYPE> v( phalf / std::sqrt( 1.0 + phalf.dot(phalf) ) );
    //   Scalar vv = v.dot( v );
    //   Vec3<MOM_TYPE> a( dp / dt ); // a is for now force, will be converted to dv/dt
    //   // convert a to dv/dt
    //   a = ( a - v * ( v.dot(a) ) ) * std::sqrt( 1.0 - vv );
    //   Scalar va = v.dot( a ); // get the real v dot a
    //   return vv / std::max( std::sqrt( a.dot(a) - va * va / vv ), 1e-6 ); // in case denominator becomes zero
    // }

    // float ParticlePusher::GetDipolarRc(const Scalar &r_sph, const Scalar &cos_th, const Scalar &phi) const {
    //   Scalar sin_th = std::sqrt( 1.0 - cos_th * cos_th );
    //   Scalar tmp1 = 1.0 + cos_th * cos_th;
    //   Scalar tmp2 = 3.0 * tmp1 - 2.0;
    //   return r_sph * tmp2 * std::sqrt(tmp2) / ( 3.0 * tmp1 * sin_th );
    // }

    T (*calc_Rc) ( const Ptc& ptc, const apt::Vec<T,Ptc::NDim>& dp, T dt ) = nullptr;
    T (*sample_E_ph) () = [](){ return 3.5;};

    std::optional<T> operator() ( const Ptc& ptc, const apt::Vec<T,Ptc::NDim>& dp, T dt,
                                  const apt::Vec<T,Ptc::NDim>& B, util::Rng<T>& rng ) override {
      T Rc = calc_Rc( ptc, dp, dt );
      T gamma = std::sqrt( (is_massive) + apt::sqabs(ptc.p()) );

      if (gamma > gamma_off && gamma > K_thr *  std::cbrt(Rc) && rng.uniform() < prob_fiducial * gamma  / Rc ) {
        // return sampled energy
        return {std::min( sample_E_ph(), gamma - 1.0) };
      }
      return {};

    }
  };

  template < class Ptc >
  struct MagneticConvert : public scat::Channel<Ptc> {
    using T = typename Ptc::vec_type::element_type;
    T B_thr = 1000;
    T mfp = 0.2;

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
    T mfp = 5.0;

    std::optional<T> operator() ( const Ptc& photon, const apt::Vec<T,Ptc::NDim>& dp, T dt,
                                  const apt::Vec<T,Ptc::NDim>& B, util::Rng<T>& rng ) override {
      // prob_mag_conv = dt / mfp_ph_ph
      return ( rng.uniform() < dt / mfp )  ? std::optional<T>(0.0) : std::nullopt;
    }
  };

}

namespace particle {
  template < typename Real, int DGrid, int DPtc, typename state_t, typename ShapeF,
             typename Real_dJ, knl::coordsys CS >
  ParticleUpdater< Real, DGrid, DPtc, state_t, ShapeF, Real_dJ, CS >
  ::ParticleUpdater( const knl::Grid< Real, DGrid >& localgrid, const util::Rng<Real>& rng, const std::optional<mpi::CartComm>& cart, const dye::Ensemble<DGrid>& ensemble )
    : _localgrid(localgrid),
      _dJ( knl::dims(localgrid), ShapeF() ),
      _rng(rng), _cart(cart), _ensemble(ensemble) {
    // TODOL
    using PtcArr = array<Real, DPtc, state_t>;
    using Ptc = typename PtcArr::particle_type;
    force::Factory<Ptc>::Register("lorentz", force::lorentz<Ptc>);
    force::Factory<Ptc>::Register("landau0", force::landau0<Ptc>);

    {
      scat::RadiationFromCharges<false,PtcArr> ep_scat;
      ep_scat.add(InZone<Ptc>{});
      ep_scat.add(CurvatureRadiate<Ptc>{});

      scat_gen<PtcArr>.Register( species::electron, ep_scat );
      scat_gen<PtcArr>.Register( species::positron, ep_scat );
    }

    {
      scat::PhotonPairProduction<PtcArr>photon_scat;
      photon_scat.add(InZone<Ptc>{});
      photon_scat.add(MagneticConvert<Ptc>{});
      photon_scat.add(TwoPhotonCollide<Ptc>{});
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

    std::vector<force::specs<Ptc>> specs {
                                          { "lorentz", static_cast<Real>(prop.charge_x) / prop.mass_x },
                                          { "landau0", 1000 },
    };

    auto update_p_opt = std::make_optional
      ( [specs=std::move(specs)]( auto& ptc, auto&&... args ) {
        for ( const auto& spec : specs ) {
          auto f = force::Factory<Ptc>::create(spec.id);
          f( ptc, std::forward<decltype(args)>(args)..., spec.param0 );
        }
      } );

    if ( prop.charge_x == 0 ) update_p_opt.reset(); // TODOL maybe also check is_massive

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
        if ( update_p_opt ) (*update_p_opt)( ptc, dt, E_itpl, B_itpl );
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
      auto& Jmesh = _dJ.integrate( *_cart );

      // rescale Jmesh back to real grid delta
      apt::foreach<0, DGrid> // NOTE it's DGrid not DField
        ( [&]( auto comp, const auto& g ) { // comp is proxy
            auto tmp = g.delta();
            for ( auto& elm : comp.data() ) elm *= tmp;
          }, Jmesh, _localgrid );

      // J = getJfromJmesh(_dJ); TODO coordinate system
    }

  }

}

#include "traits.hpp"
using namespace traits;
namespace particle {
  template class ParticleUpdater< real_t, DGrid, DPtc, ptc_state_t, ShapeF,
                                  real_dj_t, coordinate_system >;
}
