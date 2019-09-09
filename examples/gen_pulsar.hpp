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

#include "ic/multipole.hpp"

#include "bc/axissymmetric.hpp"
#include "bc/injection.hpp"

#include "dye/ensemble.hpp"

#include "pic.hpp"

// common parameters
namespace pic {
  constexpr long double PI = std::acos(-1.0l);

  inline constexpr const char* project_name = "Pulsar256";
  inline constexpr const char* datadir_prefix = "/home/hooray/Projects/Pigeon/Data/";

  inline constexpr apt::array<int,DGrid> dims = { 1, 1 };
  inline constexpr apt::array<bool,DGrid> periodic = {false,false};
  inline constexpr real_t dt = 0.003;

  constexpr mani::Grid<real_t,DGrid> supergrid
  = {{ { 0.0, std::log(30.0), 256 }, { 0.0, PI, 256 } }};
  inline constexpr int guard = 1;

  inline constexpr real_t wdt_pic = 1.0 / 30.0;
  inline constexpr real_t w_gyro_unitB = 3750; // set the impact of unit field strength on particle

  inline void set_resume_dir( std::optional<std::string>& dir ) {}
  inline constexpr int total_timesteps = 20000;
}

namespace pic {
  inline constexpr ModuleRange sort_particles_mr { true, 0, 100 };

  inline constexpr ModuleRange export_data_mr { true, 0, 200 };
  inline constexpr int pmpio_num_files = 1;
  inline constexpr int downsample_ratio = 1;

  inline constexpr ModuleRange checkpoint_mr { false, 1, 10000 };
  inline constexpr int num_checkpoint_parts = 4;
  inline constexpr int max_num_ckpts = 2;
  inline constexpr std::optional<float> checkpoint_autosave_hourly;

  inline constexpr ModuleRange dlb_mr { false, 1, 1000 };
  inline constexpr std::optional<int (*) ( int )> dlb_init_replica_deploy {}; // take in ensemble label and return the intended number of replicas in that ensemble
  inline constexpr std::size_t dlb_target_load = 100000;

  inline constexpr ModuleRange msperf_mr { true, 0, 100 };
  inline constexpr std::optional<int> msperf_max_entries {};
  inline constexpr auto msperf_qualified =
    []( const std::optional<dye::Ensemble<DGrid>>& ens_opt ) -> bool { return true; };

  inline constexpr ModuleRange vitals_mr { true, 0, 100 };

  inline constexpr int cout_ts_interval = 10;
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

namespace pic {
  template < int DGrid, typename Real >
  auto set_up_initial_conditions( const mani::Grid<Real,DGrid>& grid ) {
    ic::MagneticDipole<DGrid> dipole;
    {
      apt::tie(dipole.Ib[0], dipole.extent[0]) = gtl( {0.0, std::log(30.0)}, grid[0] );
      apt::tie(dipole.Ib[1], dipole.extent[1]) = gtl( {0.0, PI}, grid[1] );
    }

    return dipole;
  }

  constexpr real_t N_atm_floor = std::exp(1.0) * 2.0 * field::Omega * pic::w_gyro_unitB * std::pow( dt / wdt_pic, 2.0 );
  constexpr real_t N_atm_x = 2.0;
  constexpr real_t v_th = 0.3;
  constexpr real_t gravity_strength = 1.8;

  template < int DGrid, typename Real >
  auto set_up_boundary_conditions( const mani::Grid<Real,DGrid>& grid ) {
    bc::Axissymmetric<DGrid> axis;
    {
      constexpr auto AxisDir = decltype(axis)::AxisDir;
      axis.is_lower = std::abs( grid[AxisDir].lower() - 0.0 ) < grid[AxisDir].delta();
      axis.is_upper = std::abs( grid[AxisDir].upper() - std::acos(-1.0) ) < grid[AxisDir].delta();
    }

    bc::Injector<DGrid,Real> inj;
    {
      apt::tie(inj.Ib[0], inj.extent[0]) = gtl( field::ofs::indent[0] - 1, 1, supergrid[0], grid[0] );
      apt::tie(inj.Ib[1], inj.extent[1]) = gtl( {0.0, PI}, grid[1] );
      inj.v_th = v_th;
      inj.N_atm = N_atm_x * N_atm_floor;
      inj.posion = particle::species::ion;
      inj.negaon = particle::species::electron;
      inj.omega_t = field::omega_spinup;
    }

    return std::make_tuple(axis,inj);
  }
}

namespace particle {
  constexpr pic::real_t gamma_thr = 20.0;
  constexpr pic::real_t gamma_off = 15.0;
  constexpr pic::real_t emission_rate = 0.25;
  constexpr pic::real_t E_ph = 4.0;

  // NOTE called in main. This guarantees that idles also have properties set up correctly, which is currently a requirement for doing dynamic balance correct
  template < typename Real = double > // TODOL this is an ad hoc trick to prevent calling properties in routines where it is not needed
  void set_up_properties() {
    properties.insert(species::electron, {1,-1,"electron","el"});
    properties.insert(species::positron, {1,1,"positron","po"});
    properties.insert(species::ion, { 5, 1, "ion","io"});
    properties.insert(species::photon, { 0, 0, "photon","ph" });
  }

  // NOTE called in particle updater
  template < typename Real >
  void set_up() {
    using namespace pic;
    {
      constexpr auto* lorentz = force::template lorentz_exact<real_t,Specs,vParticle>;
      constexpr auto* landau0 = force::landau0<real_t,Specs,vParticle>;
      constexpr auto* gravity = force::gravity<real_t,Specs,vParticle>;
      real_t landau0_B_thr = 0.1;

      if ( properties.has(species::electron) ) {
        auto sp = species::electron;
        Force<real_t,Specs> force;
        const auto& prop = properties[sp];

        force.add( lorentz, ( pic::w_gyro_unitB * prop.charge_x ) / prop.mass_x );
        force.add( gravity, pic::gravity_strength );
        // force.add( landau0, landau0_B_thr );

        force.Register(sp);
      }
      if ( properties.has(species::positron) ) {
        auto sp = species::positron;
        Force<real_t,Specs> force;
        const auto& prop = properties[sp];

        force.add( lorentz, ( pic::w_gyro_unitB * prop.charge_x ) / prop.mass_x );
        force.add( gravity, pic::gravity_strength );
        // force.add( landau0, landau0_B_thr );

        force.Register(sp);
      }
      if ( properties.has(species::ion) ) {
        auto sp = species::ion;
        Force<real_t,Specs> force;
        const auto& prop = properties[sp];

        force.add( lorentz, ( pic::w_gyro_unitB * prop.charge_x ) / prop.mass_x );
        force.add( gravity, pic::gravity_strength );
        // force.add( landau0, landau0_B_thr );

        force.Register(sp);
      }
    }

    using Ptc_t = typename array<real_t,Specs>::particle_type;
    {
      Scat<real_t,Specs> ep_scat;

      ep_scat.eligs.push_back([](const Ptc_t& ptc){ return ptc.q()[0] < std::log(9.0); });

      scat::CurvatureRadiation<real_t,Specs>::K_thr = gamma_thr;
      scat::CurvatureRadiation<real_t,Specs>::gamma_off = gamma_off;
      scat::CurvatureRadiation<real_t,Specs>::emission_rate = emission_rate;
      scat::CurvatureRadiation<real_t,Specs>::sample_E_ph = []() noexcept -> real_t { return E_ph; };
      ep_scat.channels.push_back( scat::CurvatureRadiation<real_t,Specs>::test );

      if ( properties.has(species::photon) )
        ep_scat.impl = scat::RadiationFromCharges<false,real_t,Specs>;
      else
        ep_scat.impl = scat::RadiationFromCharges<true,real_t,Specs>;

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

#include <sstream>
#include "apt/print.hpp"

namespace pic {
  constexpr real_t classic_electron_radius () noexcept {
    real_t res = wdt_pic * wdt_pic / ( 4.0 * std::acos(-1.0l) * dt * dt);
    apt::foreach<0,DGrid>( [&res](const auto& g) { res *= g.delta(); }, supergrid );
    return res;
  }

  std::string characteristics(std::string indent) {
    std::ostringstream o;
    auto gamma_0 = std::pow(field::Omega,2.0) * w_gyro_unitB;
    o << indent << "gamma_0=" << apt::fmt("%.0f", gamma_0 ) << std::endl;
    o << indent << "w_pic dt=" << apt::fmt("%.4f", wdt_pic ) << std::endl;
    o << indent << "Ndot_GJ=" << apt::fmt("%.4e", gamma_0 / classic_electron_radius() ) << std::endl;

    o << indent << "ATM: N_atm_floor=" << apt::fmt("%.1f", N_atm_floor);
    o << ", multiplicity=" << apt::fmt("%.1f", N_atm_x);
    o << ", v_th=" << apt::fmt("%.2f", v_th) << ", g=" << apt::fmt("%.2f", gravity_strength) << std::endl;

    o << indent << "PC: gamma_thr=" << apt::fmt("%.0f", particle::gamma_thr);
    o << ", gamma_off=" << apt::fmt("%.0f", particle::gamma_off);
    o << ", gamma_ph=" << apt::fmt("%.0f", particle::E_ph);
    o << ", emission_rate=" << apt::fmt("%.2f", particle::emission_rate);
    return o.str();
  }
}

#endif
