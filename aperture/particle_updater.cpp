#include "particle_updater.hpp"
#include "dye/ensemble.hpp"

#include "particle/forces.hpp"
#include "particle/pair_producer.hpp"
#include "particle/migration.hpp"

#include "kernel/coordinate.hpp"

#include "kernel/grid.hpp"
#include "kernel/shapef.hpp"

namespace particle {
  map<Properties> properties;
}

namespace particle {
  template < typename T, typename p_t, typename dp_t >
  T calc_Rc( T dt, const apt::VecExpression<p_t>& p, const apt::VecExpression<dp_t>& dp ) noexcept {
    // TODO don't use a uniform number for Rc
    return 1.0;
  }

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
}

namespace particle {
  template < typename Real, int DGrid, int DPtc, typename state_t, typename ShapeF,
             typename Real_dJ,
             particle::PairScheme pair_scheme, knl::coordsys CS >
  ParticleUpdater< Real, DGrid, DPtc, state_t, ShapeF, Real_dJ, pair_scheme, CS >
  ::ParticleUpdater( const knl::Grid< Real, DGrid >& localgrid, const util::Rng<Real>& rng, const std::optional<mpi::CartComm>& cart, const dye::Ensemble<DGrid>& ensemble )
    : _localgrid(localgrid),
      _dJ( knl::dims(localgrid), ShapeF() ),
      _rng(rng), _cart(cart), _ensemble(ensemble) {}
}

// {
//   template < typename Ptc, typename T = typename Ptc::value_type >
//     bool is_productive_lepton( const PtcExpression<Ptc>& ptc, const T& gamma,
//                                const T& Rc, util::Rng<T>& rng ) noexcept;
//   template < typename Ptc, typename T = typename Ptc::value_type >
//     bool is_productive_photon( const PtcExpression<Ptc>& photon, const T& B2,
//                                util::Rng<T>& rng ) noexcept;

// namespace particle {
//   template < typename Ptc, typename T >
//   bool is_productive_lepton( const PtcExpression<Ptc>& ptc, const T& gamma,
//                              const T& Rc, util::Rng<T>& rng ) noexcept {
//     // TODO pane
//     // return
//     //   _pane.in_productive_zone( ptc.q ) &&
//     //   gamma > _pane.gamma_off &&
//     //   gamma > _pane.K_thr *  std::cbrt(Rc) &&
//     //   // prob_fiducial = K_curv_em_rate * dt
//     //   rng.uniform() < _pane.prob_fiducial * gamma  / Rc;
//   }
// }

// namespace particle {
//   namespace opacity {
//     // prob_mag_conv = dt / mfp_mag_conv
//     template < typename T >
//     inline bool mag_conv( const T& B2, T randnum ) noexcept {
//       // TODO pane
//       // return B2 > _pane.B_magconv * _pane.B_magconv ? ( randnum < prob_mag_conv ) : false;
//     }

//     // template < typename T >
//     // inline T f_x ( T x ) {
//     //   // distribution of x*exp(-x^2/2), which peaks at x = 1.
//     //   return std::sqrt( -2.0 * std::log(x) );
//     // }

//     // prob_mag_conv = dt / mfp_ph_ph
//     template < typename T >
//     inline bool ph_ph( T randnum ) noexcept {
//       // TODO double check this implementation, it is not equivalent because the original one has some sort of gaussian in it. Use Monte Carlo
//       // TODO prob_ph_ph not defined
//       // return randnum < prob_ph_ph;
//     }
//   }

//   namespace particle {
//     template < typename T >
//     inline T sample_E_ph(const T& gamma, const T& Rc) {
//       // TODO pane
//       // return std::min(_pane.E_ph, gamma - 1.0);
//     }

//   }

//   template < typename Ptc, typename T >
//   bool is_productive_photon( const PtcExpression<Ptc>& photon, const T& B2,
//                              util::Rng<T>& rng ) noexcept {
//     return opacity::mag_conv(B2, rng.uniform() ) || opacity::ph_ph( rng.uniform() );
//   }
// }
//   if constexpr ( IsCharged && IsRadiative ) {
//       // TODOL wrap into a factory of pair producer
//       if constexpr ( pair_scheme != PairScheme::Disabled ) {
//           auto gamma = std::sqrt( (!prop.mass_x) + apt::sqabs(ptc.p()) );
//           if ( is_productive_lepton( ptc, gamma, Rc, _rng ) ) {
//             if constexpr ( pair_scheme == PairScheme::Photon ) {
//                 produce_photons( std::back_inserter( particles[species::photon] ),
//                                  ptc, gamma, Rc );
//               } else if ( pair_scheme == PairScheme::Instant ) {
//               instant_produce_pairs( std::back_inserter( particles[species::electron] ),
//                                      std::back_inserter( particles[species::positron] ),
//                                      ptc, gamma, Rc );
//             }
//           }
//         }
//     }
//   else if ( pair_scheme == PairScheme::Photon ) {
//     // TODO check is photon
//     is_productive_photon;
//     photon_produce_pairs( std::back_inserter( particles[species::electron] ),
//                           std::back_inserter( particles[species::positron] ),
//                           ptc );
//   }
// }

namespace particle {
  template < typename Real, int DGrid, int DPtc, typename state_t, typename ShapeF,
             typename Real_dJ,
             PairScheme pair_scheme, knl::coordsys CS >
  template < bool IsCharged >
  void ParticleUpdater< Real, DGrid, DPtc, state_t, ShapeF, Real_dJ, pair_scheme, CS >
  ::update_species( array<Real,3,state_t>& sp_ptcs,
                    Real dt, Real unit_e,
                    const Properties& prop,
                    const field::Field<Real,3,DGrid>& E,
                    const field::Field<Real,3,DGrid>& B,
                    const apt::array< apt::pair<Real>, DGrid >& borders
                    ) {
    if ( sp_ptcs.size() == 0 ) return;

    using Ptc = decltype(sp_ptcs[0]);
    std::vector<force::specs<Ptc>> specs {
                                          { "lorentz", static_cast<Real>(prop.charge_x) / prop.mass_x },
                                          { "landau0", 1000 },
    };

    auto update_p =
      [specs=std::move(specs)]( auto&&... args ) {
        for ( const auto& spec : specs ) {
          auto f = force::Factory<Ptc>::create(spec.id);
          f( std::forward<decltype(args)>(args)..., spec.param0 );
        }
      };

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
            return knl::coord<CS>::geodesic_move( ptc.q(), ptc.p(), dt / gamma );
          } else {
          apt::Vec<Real, DPtc> dq; // for RVO
          dq = knl::coord<CS>::geodesic_move( ptc.q(), (ptc.p() /= gamma), dt );
          ptc.p() *= gamma;
          return dq;
        }
      };

    auto abs2std =
      [&grid=_localgrid]( const auto& qabs ) {
        apt::array<Real,DGrid> q_std;
        apt::foreach<0,DGrid>
          ( [](auto& q, auto q_abs, const auto& g ) noexcept {
              q = ( q_abs - g.lower() ) / g.delta();
            }, q_std, qabs, grid );
        return q_std;
      };

    for ( auto ptc : sp_ptcs ) { // TODOL sematics, check ptc is proxy
      if( ptc.is(flag::empty) ) continue;

      // TODOL should also check IsMassive
      Real Rc{};
      if constexpr ( IsCharged ) {
        auto q_std = abs2std( ptc.q() );
        auto E_itpl = field::interpolate( E, q_std, shapef );
        auto B_itpl = field::interpolate( B, q_std, shapef );

        apt::Vec<Real,DPtc> dp = ptc.p();
        update_p( ptc, dt, E_itpl, B_itpl );
        dp -= ptc.p();
        // Rc = calc_Rc<Real>( dt, ptc.p(), std::move(dp) );
      }

      // TODO
      // pair_produce( std::back_inserter(sp_ptcs), ptc, Rc ); // use sp_ptcs temporarily to hold newly created particles if any

      // NOTE q is updated, starting from here, particles may be in the guard cells.
      auto&& dq = update_q( ptc, dt );
      {
        auto abs2std_dJ =
          [&grid=_localgrid]( apt::Vec<Real,DPtc> q1, auto dq ) {
            apt::foreach<0,DGrid> // NOTE this is DGrid, not DPtc.
              ( [](auto& q1, auto& q0, const auto& g){
                  q1 = ( q1 - g.lower() ) / g.delta();
                  q0 = q1 - q0 / g.delta();
                }, q1, dq, grid );
            return std::make_tuple(dq, q1); // NOTE the result is of DPtc
          };
        auto&&[ q0_std, q1_std ] = abs2std_dJ( ptc.q(), std::move(dq) );
        // TODOL pusher handle boundary condition. Is it needed?
        if constexpr ( IsCharged )
                       _dJ.deposit( charge_over_dt, std::move(q0_std),
                                    std::move(q1_std) ); // TODOL check the 2nd argument
      }

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

namespace particle {
  template < typename Real, int DGrid, int DPtc, typename state_t,
             typename ShapeF,
             typename Real_dJ,
             PairScheme pair_scheme, knl::coordsys CS
             >
  void ParticleUpdater< Real, DGrid, DPtc, state_t, ShapeF, Real_dJ, pair_scheme, CS >
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
      const Properties& prop = properties.at(sp);
      auto[ m_x, q_x, is_rad ] = prop;
      // TODO
      if ( q_x )
        update_species<true>( ptcs, dt, unit_e, prop, E, B, borders );
      else
        update_species<false>( ptcs, dt, unit_e, prop, E, B, borders );
      // else if ( q_x && !is_rad )
      //   update_species<true,false>( ptcs, particles, dt, unit_e, E, B, borders );
      // else
      //   update_species<false,false>( ptcs, particles, dt, unit_e, E, B, borders );
      // TODO put particles where they belong
    }

    migrate( _migrators, _ensemble.inter, timestep );
    for ( auto&& ptc : _migrators ) {
      particles[ptc.template get<species>()].push_back( std::move(ptc) ); // TODOL check this
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
                                  real_dj_t, pair_produce_scheme, coordinate_system>;
}
