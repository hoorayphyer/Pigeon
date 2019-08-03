#ifndef _GEN_HPP_
#define _GEN_HPP_

#include "apt/numeric.hpp"
#include "apt/index.hpp"

#include "pic/module_range.hpp"
#include "pic/range_conversion.hpp"
#include "pic/forces/gravity.hpp"
#include "pic/forces/landau0.hpp"

#include "manifold/grid.hpp"

#include "field/field.hpp"
#include "field/params.hpp"

#include "particle/properties.hpp"
#include "particle/map.hpp"
#include "particle/array.hpp"
#include "particle/particle.hpp"

#include "bc/axissymmetric.hpp"

#include "dye/ensemble.hpp"

#include "pic.hpp"

// common parameters
namespace pic {
  constexpr long double PI = std::acos(-1.0l);

  inline constexpr const char* project_name = "Pulsar256";
  inline constexpr const char* datadir_prefix = "/home/hooray/Projects/Pigeon/Data/";

  inline constexpr apt::array<int,DGrid> dims = { 1, 1 };
  inline constexpr apt::array<bool,DGrid> periodic = {false,false};
  inline constexpr int total_timesteps = 20000;
  inline constexpr real_t dt = 0.003;

  constexpr mani::Grid<real_t,DGrid> supergrid
  = {{ { 0.0, std::log(30.0), 256 }, { 0.0, PI, 256 } }};
  inline constexpr int guard = 1;

  inline constexpr real_t wdt_pic = 1.0 / 30.0;
  inline constexpr real_t w_gyro_unitB = 3750; // set the impact of unit field strength on particle
}

namespace pic {
  inline constexpr ModuleRange sort_particles_mr { true, 0, 100 };

  inline constexpr ModuleRange export_data_mr { true, 0, 200 };
  inline constexpr int pmpio_num_files = 1;
  inline constexpr int downsample_ratio = 1;

  inline constexpr ModuleRange checkpoint_mr { false, 1, 10000 };
  inline constexpr int num_checkpoint_parts = 4;
  inline constexpr std::optional<float> checkpoint_hourly;

  inline constexpr ModuleRange dlb_mr { false, 1, 1000 };
  inline constexpr std::optional<int (*) ( int )> dlb_init_replica_deploy {}; // take in ensemble label and return the intended number of replicas in that ensemble
  inline constexpr std::size_t dlb_target_load = 100000;

  inline constexpr ModuleRange msperf_mr { true, 0, 100 };
  inline constexpr std::optional<int> msperf_max_entries {};
  inline constexpr auto msperf_qualified =
    []( const std::optional<dye::Ensemble<DGrid>>& ens_opt ) -> bool { return true; };

  inline constexpr ModuleRange stats_mr { true, 0, 100 };

  inline constexpr int cout_ts_interval = 100;
}

namespace field {
  using pic::real_t;
  // TODOL these will become free parameters
  inline constexpr real_t Omega = 1.0 / 6.0;

  constexpr real_t omega_spinup ( real_t time ) noexcept {
    return std::min<real_t>( time / 4.0, 1.0 ) * Omega;
  }

  template < typename Real >
  void set_up() {
    ofs::magnetic_pole = 2; // 1 for mono-, 2 for di-
    ofs::indent = { 5, 43, pic::guard, pic::guard };
    ofs::damping_rate = 10.0;
  }
}

namespace particle {
  // NOTE called in main. This guarantees that idles also have properties set up correctly, which is currently a requirement for doing dynamic balance correct
  template < typename Real = double > // TODOL this is an ad hoc trick to prevent calling properties in routines where it is not needed
  void set_up_properties() {
    properties[species::electron] = {1,-1,"electron"};
    properties[species::positron] = {1,1,"positron"};
    properties[species::ion] = { 5, 1, "ion"};
    properties[species::photon] = { 0, 0, "photon" };
  }

  // NOTE called in particle updater
  template < typename Real >
  void set_up() {
    using namespace pic;
    {
      constexpr auto* lorentz = force::template lorentz<real_t,Specs,vParticle>;
      constexpr auto* landau0 = force::landau0<real_t,Specs,vParticle>;
      constexpr auto* gravity = force::gravity<real_t,Specs,vParticle>;
      real_t landau0_B_thr = 0.1;
      real_t gravity_strength = 0.5;

      if ( properties.has(species::electron) ) {
        auto sp = species::electron;
        Force<real_t,Specs> force;
        const auto& prop = properties.at(sp);

        force.add( lorentz, ( pic::w_gyro_unitB * prop.charge_x ) / prop.mass_x );
        force.add( gravity, gravity_strength );
  //      force.add( landau0, landau0_B_thr );

        force.Register(sp);
      }
      if ( properties.has(species::positron) ) {
        auto sp = species::positron;
        Force<real_t,Specs> force;
        const auto& prop = properties.at(sp);

        force.add( lorentz, ( pic::w_gyro_unitB * prop.charge_x ) / prop.mass_x );
        force.add( gravity, gravity_strength );
    //    force.add( landau0, landau0_B_thr );

        force.Register(sp);
      }
      if ( properties.has(species::ion) ) {
        auto sp = species::ion;
        Force<real_t,Specs> force;
        const auto& prop = properties.at(sp);

        force.add( lorentz, ( pic::w_gyro_unitB * prop.charge_x ) / prop.mass_x );
        force.add( gravity, gravity_strength );
        // force.add( landau0, landau0_B_thr );

        force.Register(sp);
      }
    }

    using Ptc_t = typename array<real_t,Specs>::particle_type;
    {
      Scat<real_t,Specs> ep_scat;

      ep_scat.eligs.push_back([](const Ptc_t& ptc){ return ptc.q()[0] < std::log(9.0); });

      scat::CurvatureRadiation<real_t,Specs>::K_thr = 10.0;
      scat::CurvatureRadiation<real_t,Specs>::gamma_off = 5.0;
      scat::CurvatureRadiation<real_t,Specs>::emission_rate = 0.25;
      scat::CurvatureRadiation<real_t,Specs>::sample_E_ph = []() noexcept -> real_t { return 3.0; };
      ep_scat.channels.push_back( scat::CurvatureRadiation<real_t,Specs>::test );

      ep_scat.impl = scat::RadiationFromCharges<false,real_t,Specs>;

      if ( properties.has(species::electron) && properties.has(species::positron) ) {
        ep_scat.Register( species::electron );
        ep_scat.Register( species::positron );
      }
    }

    if ( properties.has(species::photon) ) {
      Scat<real_t,Specs> photon_scat;
      // Photons are free to roam across all domain. They may produce pairs outside light cylinder
      photon_scat.eligs.push_back([](const Ptc_t& ptc) { return true; });
      scat::MagneticConvert<real_t,Specs>::B_thr = 0.1;
      scat::MagneticConvert<real_t,Specs>::mfp = 0.2;
      photon_scat.channels.push_back( scat::MagneticConvert<real_t,Specs>::test );

      scat::TwoPhotonCollide<real_t,Specs>::mfp = 5.0;;
      photon_scat.channels.push_back( scat::TwoPhotonCollide<real_t,Specs>::test );

      photon_scat.impl = scat::PhotonPairProduction<real_t,Specs>;

      photon_scat.Register( species::photon );
    }
  }

}

namespace pic {
  template < int DGrid,
             typename Real,
             template < typename > class Specs,
             typename RealJ >
  struct InitialCondition {
  private:
    const mani::Grid<Real,DGrid>& _grid;
    field::Field<Real, 3, DGrid>& _Bfield;

    apt::Index<DGrid> _Ib;
    apt::Index<DGrid> _extent;

    Real B_r ( Real logr, Real theta ) noexcept {
      return 2.0 * std::cos(theta) * std::exp(-3.0 * logr);
    };

    Real B_th ( Real logr, Real theta ) noexcept {
      return std::sin(theta) * std::exp(-3.0 * logr);
    };

  public:
    InitialCondition ( const mani::Grid<Real,DGrid>& localgrid,
                       const field::Field<Real, 3, DGrid>& Efield,
                       field::Field<Real, 3, DGrid>& Bfield,
                       const field::Field<RealJ, 3, DGrid>& Jfield, // J is Jmesh on a replica
                       const particle::map<particle::array<Real,Specs>>& particles )
      : _grid(localgrid), _Bfield(Bfield) {
      apt::tie(_Ib[0], _extent[0]) = gtl( {0.0, std::log(30.0)}, localgrid[0] );
      apt::tie(_Ib[1], _extent[1]) = gtl( {0.0, PI}, localgrid[1] );
    }

    void operator() () {
      for ( auto I : apt::Block(_extent) ) {
        I += _Ib;
        _Bfield[0](I) = B_r( _grid[0].absc(I[0], _Bfield[0].offset()[0]), _grid[1].absc(I[1], _Bfield[0].offset()[1]) );
        _Bfield[1](I) = B_th( _grid[0].absc(I[0], _Bfield[1].offset()[0]), _grid[1].absc(I[1], _Bfield[1].offset()[1]) );
      }
    }

    constexpr int initial_timestep() noexcept { return 0; }
  };

  // as a boundary condition for particles
  template < int DGrid,
             typename Real,
             template < typename > class Specs,
             typename RealJ >
  struct Injector {
  private:
    const mani::Grid<Real,DGrid>& _grid;
    const field::Field<Real, 3, DGrid>& _Efield;
    const field::Field<Real, 3, DGrid>& _Bfield;
    const field::Field<RealJ, 3, DGrid>& _Jfield; // J is Jmesh on a replica
    particle::map<particle::array<Real,Specs>>& _particles;

    apt::Index<DGrid> _Ib;
    apt::Index<DGrid> _extent;

  public:
    Injector ( const mani::Grid<Real,DGrid>& localgrid,
               const field::Field<Real, 3, DGrid>& Efield,
               const field::Field<Real, 3, DGrid>& Bfield,
               const field::Field<RealJ, 3, DGrid>& Jfield, // J is Jmesh on a replica
               particle::map<particle::array<Real,Specs>>& particles )
      : _grid(localgrid), _Efield(Efield), _Bfield(Bfield), _Jfield(Jfield), _particles(particles) {
      apt::tie(_Ib[0], _extent[0]) = gtl( field::ofs::indent[0] - 1, 1, supergrid[0], localgrid[0] );
      apt::tie(_Ib[1], _extent[1]) = gtl( {0.0, PI}, localgrid[1] );
    }

    void operator() ( int timestep, Real dt, util::Rng<Real>& rng, particle::birthplace birth ) {
      using namespace particle;

      constexpr Real v_th = 0.3;
      constexpr Real j_reg_x = 0.0;
      constexpr Real N_dot = 500.0;

      constexpr auto posion = species::ion;
      constexpr auto negaon = species::electron;

      Real omega = field::omega_spinup( timestep * dt );

      int charge_x = particle::properties.at(posion).charge_x;

      // the idea is that (timestep + 1) * N_inj * profile represents the accummulated number of injected pairs through the specified timestep. NOTE the timestep is shifted by one to reflect the actual times the inject is called, including this time.
      // NOTE _Jfield and Jmesh also differs in their underlying mesh guard cells
      auto profile_inj =
        [N_inj=N_dot*dt,timestep] ( Real theta ) noexcept {
          Real inj_num_base = N_inj * 0.5 * std::abs( std::sin( 2 * theta ) );
          return static_cast<int>( (timestep + 1) * inj_num_base ) - static_cast<int>( timestep * inj_num_base );
        };

      auto j_reg_inj =
        [&grid=_grid, charge_x, j_reg_x]( Real J ) noexcept {
          static auto factor = j_reg_x * grid[0].delta() * grid[1].delta() / charge_x;
          return J * factor;
        };

      auto itr_po = std::back_inserter(_particles[posion]);
      auto itr_ne = std::back_inserter(_particles[negaon]);

      for ( auto I : apt::Block(_extent) ) {
        I += _Ib;
        apt::Vec<Real, Specs<Real>::Dim> q{};
        for ( int i = 0; i < DGrid; ++i )
          q[i] = _grid[i].absc(I[i], 0.5);

        apt::Vec<Real,3> nB, J;
        // TODOL use interpolated value
        for ( int i = 0; i < 3; ++i ) {
          nB[i] = _Bfield[i](I);
          J[i] = _Jfield[i](I);
        }

        int num = std::max<Real>( profile_inj(q[1]), j_reg_inj(apt::abs(J)) );

        // find n_B
        nB /= apt::abs(nB);

        apt::Vec<Real, Specs<Real>::Dim> p{};
        p[2] = omega * std::exp(q[0]) * std::sin(q[1]); // corotating

        for ( int n = 0; n < num; ++n ) {
          auto q_ptc = q;
          for ( int i = 0; i < DGrid; ++i )
            q_ptc[i] += _grid[i].delta() * rng.uniform(-0.5, 0.5);
          auto p_ptc = p;
          p_ptc += nB * rng.gaussian( 0.0, v_th );
          *(itr_ne++) = Particle<Real,Specs>( q_ptc, p_ptc, negaon, birth );
          *(itr_po++) = Particle<Real,Specs>( std::move(q_ptc), std::move(p_ptc), posion, birth );
        }
      }
    }
  };
}

#endif
