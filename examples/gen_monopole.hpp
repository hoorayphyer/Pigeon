#ifndef _GEN_HPP_
#define _GEN_HPP_

#include "apt/numeric.hpp"
#include "apt/index.hpp"

#include "manifold/grid.hpp"

#include "field/field.hpp"

#include "particle/properties.hpp"
#include "particle/map.hpp"
#include "particle/array.hpp"
#include "particle/forces.hpp"
#include "particle/scattering.hpp"
#include "particle/particle.hpp"

#include "bc/fold_back_J.hpp"
#include "bc/axissymmetric.hpp"

#include "pic.hpp"

namespace particle {
  template < typename Real, template < typename > class Specs >
  ScatGen<Real, Specs> scat_gen;

  template < typename Real, template < typename > class Specs, template < typename, template < typename > class > class Ptc_t >
  ForceGen<Real, Specs,Ptc_t> force_gen;
}

// common parameters
namespace pic {
  constexpr long double PI = std::acos(-1.0l);

  inline constexpr const char* project_name = "Monopole";
  inline constexpr const char* datadir_prefix = "../Data/";

  inline constexpr apt::array<int,DGrid> dims = { 1, 1 };
  inline constexpr apt::array<bool,DGrid> periodic = {false,false};
  inline constexpr int total_timesteps = 2000;
  inline constexpr real_t dt = 0.01;

  constexpr mani::Grid<real_t,DGrid> supergrid
  = {{ { 0.0, std::log(10.0), 64 }, { 0.0, PI, 64 } }};
  inline constexpr int guard = 1;

  inline constexpr int Np = 5;
  inline constexpr real_t epsilon = 0.25;

  constexpr real_t classic_electron_radius () noexcept {
    real_t res = epsilon * epsilon / ( 4 * PI* Np * dt * dt);
    apt::foreach<0,DGrid>( [&res](const auto& g) { res *= g.delta(); }, supergrid );
    return res;
  }
}

// TODOL free parameters for specific lab
namespace pic {
  inline constexpr real_t mu0 = 15000.0;
  inline constexpr real_t omega_max = 1.0 / 6.0;

  inline constexpr int spinup_duration = 2.0;

  constexpr real_t omega_spinup ( real_t time ) noexcept {
    return std::min( time / pic::spinup_duration, 1.0 ) * pic::omega_max;
  }

  namespace ofs {
    inline constexpr int magnetic_pole = 1; // 1 for mono-, 2 for di-
    inline constexpr int indent[4] = { 5, 20, guard, guard };
    inline constexpr real_t damping_rate = 0.01;
  }
}

namespace pic {
  inline constexpr int pmpio_num_files = 1;
  inline constexpr int data_export_init_ts = 0;

  inline constexpr int checkpoint_init_ts = 0;
  inline constexpr int num_ckpt_parts = 4;

  inline constexpr int dlb_init_ts = 0;
  inline constexpr std::size_t dlb_target_load = 100000;
}

namespace pic :: interval {
  inline constexpr int data_export = 20;
  inline constexpr int checkpoint = 10000;
  inline constexpr int dlb = 1000;
}

// TODOL all the stuff under this {} are meant to be user-specified. Here the pulsar in LogSpherical is used
namespace particle {
  template < typename T, template < typename > class Specs >
  using Pat = typename array< T, Specs >::particle_type; // type of particle from array

  // LogSpherical
  namespace force {
    template < typename T, template < typename > class Specs, template < typename, template < typename > class > class Ptc_t >
    void gravity( Ptc_t<T,Specs>& ptc, T dt, const apt::Vec<T,Specs<T>::Dim>& , const apt::Vec<T,Specs<T>::Dim>&, T g  ) noexcept {
      ptc.p()[0] -= g * std::exp( - 2 * ptc.q()[0] ) * dt;
    };
  }
}

namespace particle {
  template < typename Real >
  void set_up() {
    using namespace pic;
    {
      // TODO initialize properties at a different place
      // TODO move out
      properties[species::electron] = {1,-1,"electron"};
      properties[species::positron] = {1,1,"positron"};
      properties[species::ion] = { 5, 1, "ion"};
    }

    {
      constexpr Real gravity_strength = 1.8;
      using Force = force::Force<real_t,Specs,vParticle>;
      constexpr auto& fgen = force_gen<real_t,Specs,vParticle>;
      constexpr auto* lorentz = force::template lorentz<real_t,Specs,vParticle>;
      constexpr auto* gravity = force::gravity<real_t,Specs,vParticle>;
      {
        auto sp = species::electron;
        Force force;
        const auto& prop = properties.at(sp);

        force.add( lorentz, static_cast<Real>(prop.charge_x) / prop.mass_x );
        force.add( gravity, gravity_strength );

        fgen.Register( sp, force );
      }
      {
        auto sp = species::positron;
        Force force;
        const auto& prop = properties.at(sp);

        force.add( lorentz, static_cast<Real>(prop.charge_x) / prop.mass_x );
        force.add( gravity, gravity_strength );

        fgen.Register( sp, force );
      }
      {
        auto sp = species::ion;
        Force force;
        const auto& prop = properties.at(sp);

        force.add( lorentz, static_cast<Real>(prop.charge_x) / prop.mass_x );
        force.add( gravity, gravity_strength );

        fgen.Register( sp, force );
      }
    }

  }
}

namespace pic {
  // TODOL check gtl boundary on extent.
  // NOTE range is assumed to be [,), and range[1] >= range[0]
  template < typename Real >
  apt::pair< int > gtl ( const apt::pair<Real>& range,
                          const mani::Grid1D<Real>& localgrid ) noexcept {
    auto to_grid_index =
      [&localgrid] ( auto absc ) noexcept {
        return ( absc - localgrid.lower() ) / localgrid.delta();
      };
    int lb = to_grid_index( range[LFT] );
    int ub = to_grid_index( range[RGT] );

    if ( ub <= 0 || lb >= localgrid.dim() ) {
      // range not applicable on current local patch
      return {0,0};
    } else {
      lb = std::max<int>( 0, lb );
      ub = std::min<int>( ub, localgrid.dim() );
      return { lb, ub - lb };
    }
  }

  template < typename Real >
  apt::pair< int > gtl ( int Ib_global, int extent_global,
                         const mani::Grid1D<Real>& supergrid,
                         const mani::Grid1D<Real>& localgrid ) noexcept {
    // lb is the Ib_global with respect to the localgrid
    int lb = Ib_global - static_cast<int>( ( localgrid.lower() - supergrid.lower() ) / localgrid.delta() + 0.5 );
    int ub = lb + extent_global;

    if ( ub <= 0 || lb >= localgrid.dim() ) {
      // range not applicable on current local patch
      return {0,0};
    } else {
      lb = std::max<int>( 0, lb );
      ub = std::min<int>( ub, localgrid.dim() );
      return { lb, ub - lb };
    }
  }

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

    const Real _mu0 = pic::mu0;

    Real B_r_over_mu0 ( Real logr ) noexcept {
      return std::exp(-2.0 * logr);
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
        _Bfield[0](I) = _mu0 * B_r_over_mu0( _grid[0].absc(I[0], _Bfield[0].offset()[0]) );
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
      apt::tie(_Ib[0], _extent[0]) = gtl( pic::ofs::indent[0] - 1, 1, supergrid[0], localgrid[0] );
      apt::tie(_Ib[1], _extent[1]) = gtl( {0.0, PI}, localgrid[1] );
    }

    void operator() ( int timestep, Real dt, util::Rng<Real>& rng ) {
      using namespace particle;

      constexpr Real v_th = 0.3;
      constexpr Real Ninj = 1;

      constexpr auto posion = species::positron;
      constexpr auto negaon = species::electron;

      Real omega = pic::omega_spinup( timestep * dt );

      int charge_x = particle::properties.at(posion).charge_x;

      // the idea is that (timestep + 1) * N_inj * profile represents the accummulated number of injected pairs through the specified timestep. NOTE the timestep is shifted by one to reflect the actual times the inject is called, including this time.
      // NOTE _Jfield and Jmesh also differs in their underlying mesh guard cells
      auto profile_inj =
        [&Ninj,&timestep] ( Real theta ) noexcept {
          Real inj_num_base = Ninj * std::abs( std::sin(theta) );
          return static_cast<int>( (timestep + 1) * inj_num_base ) - static_cast<int>( timestep * inj_num_base );
        };

      auto itr_po = std::back_inserter(_particles[posion]);
      auto itr_ne = std::back_inserter(_particles[negaon]);

      for ( auto I : apt::Block(_extent) ) {
        I += _Ib;
        apt::Vec<Real, Specs<Real>::Dim> q{};
        for ( int i = 0; i < DGrid; ++i )
          q[i] = _grid[i].absc(I[i]);

        apt::Vec<Real,3> nB, J;
        for ( int i = 0; i < 3; ++i ) {
          nB[i] = _Bfield[i](I);
          J[i] = _Jfield[i](I);
        }

        int num = profile_inj(q[1]);

        // find n_B
        nB /= apt::abs(nB);

        apt::Vec<Real, Specs<Real>::Dim> p{};
        p[2] = omega * std::exp(q[0]) * std::sin(q[1]); // corotating

        for ( int n = 0; n < num; ++n ) {
          auto q_ptc = q;
          for ( int i = 0; i < DGrid; ++i )
            q_ptc[i] += _grid[i].delta() * rng.uniform();
          auto p_ptc = p;
          p_ptc += nB * rng.gaussian( 0.0, v_th );
          *(itr_ne++) = Particle<Real,Specs>( q_ptc, p_ptc, negaon );
          *(itr_po++) = Particle<Real,Specs>( std::move(q_ptc), std::move(p_ptc), posion );
        }
      }
    }
  };
}

#endif
