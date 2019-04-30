#ifndef _GEN_HPP_
#define _GEN_HPP_

#include "apt/numeric.hpp"
#include "apt/index.hpp"

#include "kernel/grid.hpp"

#include "field/field.hpp"

#include "particle/map.hpp"
#include "particle/array.hpp"
#include "particle/forces.hpp"
#include "particle/scattering.hpp"
#include "particle/particle.hpp"

#include "pic.hpp"

namespace particle {
  struct Properties {
    unsigned int mass_x = 0; // in terms of unit mass
    int charge_x = 0; // in terms of unit charge
    const std::string name  = "";
  };
}

namespace particle {
  extern map<Properties> properties;

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
  inline constexpr int total_timesteps = 1000;
  inline constexpr real_t dt = 0.005;

  constexpr knl::Grid<real_t,DGrid> supergrid
  = {{ { 0.0, std::log(30.0), 128 }, { 0.0, PI, 128 } }};
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

  inline constexpr int spinup_duration = 10.0;

  constexpr real_t omega_spinup ( real_t time ) noexcept {
    return std::min( time / pic::spinup_duration, 1.0 ) * pic::omega_max;
  }

  namespace ofs {
    inline constexpr int magnetic_pole = 1; // 1 for mono-, 2 for di-
    inline constexpr int indent[4] = { 5, 20, guard, guard };
    inline constexpr real_t damping_rate = 0.01;
  }
}

namespace pic :: interval {
  inline constexpr int data_export = 20;
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
  template < typename Real >
  apt::pair< Real > gtl ( const apt::pair<Real>& range,
                          const knl::Grid1D<Real>& localgrid ) noexcept {
    auto to_grid_index =
      [&localgrid] ( auto absc ) noexcept {
        return ( absc - localgrid.lower() ) / localgrid.delta();
      };
    Real Ib = std::max<int>( 0, to_grid_index( range[LFT] ) );
    Real extent = std::min<int>( to_grid_index( range[RGT] ), localgrid.dim() ) - Ib;
    return { Ib, extent };
  }

  template < int DGrid,
             typename Real,
             template < typename > class Specs,
             typename RealJ >
  struct InitialCondition {
  private:
    const knl::Grid<Real,DGrid>& _grid;
    field::Field<Real, 3, DGrid>& _Bfield;

    apt::Index<DGrid> _Ib;
    apt::Index<DGrid> _extent;

    const Real _mu0 = pic::mu0;

    Real B_r_over_mu0 ( Real logr ) noexcept {
      return std::exp(-2.0 * logr);
    };

  public:
    InitialCondition ( const knl::Grid<Real,DGrid>& localgrid,
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

  // template < int DGrid,
  //            typename Real,
  //            template < typename > class Specs,
  //            typename RealJ >
  // struct FieldBC_Rotating_Conductor {
  // private:
  //   const knl::Grid<Real,DGrid>& _grid;
  //   field::Field<Real, 3, DGrid>& _Efield;
  //   field::Field<Real, 3, DGrid>& _Bfield;

  //   apt::Index<DGrid> _Ib;
  //   apt::Index<DGrid> _extent;

  //   const Real _mu0 = pic::mu0;

  //   Real B_r_over_mu0 ( Real logr, Real theta ) noexcept {
  //     return 2.0 * std::cos(theta) * std::exp(-3.0 * logr);
  //   };

  //   Real B_th_over_mu0 ( Real logr, Real theta ) noexcept {
  //     return std::sin(theta) * std::exp(-3.0 * logr);
  //   };

  //   // find out E by E = -( Omega x r ) x B
  //   Real E_r_over_mu0omega ( Real logr, Real theta ) noexcept {
  //     auto sin_t = std::sin(theta);
  //     return std::exp( -2.0 * logr ) * sin_t * sin_t;
  //   };

  //   Real E_th_over_mu0omega ( Real logr, Real theta ) noexcept {
  //     return - std::exp( -2.0 * logr ) * std::sin( 2.0 * theta );
  //   };

  // public:
  //   FieldBC_Rotating_Conductor ( const knl::Grid<Real,DGrid>& localgrid,
  //                                field::Field<Real, 3, DGrid>& Efield,
  //                                field::Field<Real, 3, DGrid>& Bfield,
  //                                const field::Field<RealJ, 3, DGrid>& Jfield, // J is Jmesh on a replica
  //                                const particle::map<particle::array<Real,Specs>>& particles )
  //     : _grid(localgrid), _Efield(Efield), _Bfield(Bfield) {
  //     _Ib[0] = 0;
  //     _extent[0] = pic::ofs::indent[0];
  //     apt::tie(_Ib[1], _extent[1]) = gtl( {0.0, PI}, localgrid[1] );
  //   }

  //   void operator() ( int timestep, Real dt ) {
  //     const auto omega = pic::omega_spinup( timestep * dt );
  //     for ( auto I : apt::Block(_extent) ) {
  //       I += _Ib;
  //       // TODO deal with discontinuous and continuous variables
  //       // TODO piecing together different patches
  //       _Bfield[0](I) = _mu0 * B_r_over_mu0( _grid[0].absc(I[0], _Bfield[0].offset()[0]), _grid[1].absc(I[1], _Bfield[0].offset()[1]) );
  //       _Bfield[1](I) = _mu0 * B_th_over_mu0( _grid[0].absc(I[0], _Bfield[1].offset()[0]), _grid[1].absc(I[1], _Bfield[1].offset()[1]) );
  //       _Bfield[2](I) = 0.0;

  //       _Efield[0](I) = _mu0 * omega * E_r_over_mu0omega( _grid[0].absc(I[0], _Efield[0].offset()[0]), _grid[1].absc(I[1], _Efield[0].offset()[1]) );
  //       _Efield[1](I) = _mu0 * omega * E_th_over_mu0omega( _grid[0].absc(I[0], _Efield[1].offset()[0]), _grid[1].absc(I[1], _Efield[1].offset()[1]) );
  //       _Efield[2](I) = 0.0;
  //     }
  //   }
  // };

  // template < bool IsLower, int DGrid,
  //            typename Real,
  //            template < typename > class Specs,
  //            typename RealJ >
  // struct FieldBC_Axis {
  // private:
  //   const knl::Grid<Real,DGrid>& _grid;
  //   field::Field<Real, 3, DGrid>& _Efield;
  //   field::Field<Real, 3, DGrid>& _Bfield;

  //   bool _is_at_axis = false;

  // public:
  //   FieldBC_Axis ( const knl::Grid<Real,DGrid>& localgrid,
  //                       field::Field<Real, 3, DGrid>& Efield,
  //                       field::Field<Real, 3, DGrid>& Bfield,
  //                       const field::Field<RealJ, 3, DGrid>& Jfield,
  //                       const particle::map<particle::array<Real,Specs>>& particles )
  //     : _grid(localgrid), _Efield(Efield), _Bfield(Bfield) {
  //     if constexpr ( IsLower )
  //                    _is_at_axis = std::abs( localgrid[1].lower() - 0.0 ) < localgrid[1].delta();
  //     else
  //       _is_at_axis = std::abs( localgrid[1].upper() - PI ) < localgrid[1].delta();
  //   }

  //   void operator() () {
  //     // TODO We don't need to do anything to the guard cells right?
  //     // TODO NO! Guard cells values are needed when doing interpolating E and B
  //     if ( !_is_at_axis ) return;
  //     // E_theta, B_r, B_phi are on the axis. All but B_r should be set to zero
  //     const auto& mesh = _Efield.mesh();
  //     int n = IsLower ? 0 : _grid[1].dim();
  //     for ( const auto& trI : mesh.project(1, {}, mesh.extent() ) ) {
  //       _Efield[1][trI | n] = 0.0;
  //       _Bfield[2][trI | n] = 0.0;
  //     }
  //   }
  // };

  template < bool IsLower, int DGrid,
             typename Real,
             template < typename > class Specs,
             typename RealJ >
  struct FieldBC_FoldBackJ {
  private:
    static_assert(DGrid==2);
    const knl::Grid<Real,DGrid>& _grid;
    field::Field<RealJ, 3, DGrid>& _Jmesh;

    bool _is_at_axis = false;

  public:
    FieldBC_FoldBackJ ( const knl::Grid<Real,DGrid>& localgrid,
                        const field::Field<Real, 3, DGrid>& Efield,
                        const field::Field<Real, 3, DGrid>& Bfield,
                        field::Field<RealJ, 3, DGrid>& Jfield,
                        const particle::map<particle::array<Real,Specs>>& particles )
      : _grid(localgrid), _Jmesh(Jfield) {
      if constexpr ( IsLower )
                     _is_at_axis = std::abs( localgrid[1].lower() - 0.0 ) < localgrid[1].delta();
      else
        _is_at_axis = std::abs( localgrid[1].upper() - PI ) < localgrid[1].delta();
    }

    void operator() () {
      if ( !_is_at_axis ) return;

      const auto& mesh = _Jmesh.mesh();
      const int guard = mesh.guard();

      if constexpr ( IsLower ) {
          for ( const auto& trI : mesh.project(1, {}, mesh.extent() ) ) {
            for ( int n = 0; n < guard; ++n ) {
              _Jmesh[0][ trI | n ] += _Jmesh[0][ trI | -1 - n ];
              _Jmesh[2][ trI | n ] -= _Jmesh[0][ trI | -1 - n ];
            }

            _Jmesh[1][trI | 0] = 0.0;
            // TODO check the negative sign here
            for ( int n = 0; n < guard; ++n )
              _Jmesh[1][ trI | 1 + n ] -= _Jmesh[1][ trI | -1 - n ];
          }
        } else {

        const int dim = mesh.bulk_dim(1);
        for ( const auto& trI : mesh.project(1, {}, mesh.extent() ) ) {
          for ( int n = 0; n < guard; ++n ) {
            _Jmesh[0][ trI | dim - 1 - n ] += _Jmesh[0][ trI | dim + n ];
            _Jmesh[2][ trI | dim - 1 - n ] -= _Jmesh[0][ trI | dim + n ];
          }

          _Jmesh[1][trI | dim] = 0.0;
          for ( int n = 0; n < guard - 1; ++n ) // NOTE -1 is because of the favoring-lower convention
            _Jmesh[1][ trI | dim-1-n ] -= _Jmesh[1][ trI | dim+1+n ];

        }
      }
    }
  };

  // as a boundary condition for particles
  template < int DGrid,
             typename Real,
             template < typename > class Specs,
             typename RealJ >
  struct Injector {
  private:
    const knl::Grid<Real,DGrid>& _grid;
    const field::Field<Real, 3, DGrid>& _Efield;
    const field::Field<Real, 3, DGrid>& _Bfield;
    const field::Field<RealJ, 3, DGrid>& _Jfield; // J is Jmesh on a replica
    particle::map<particle::array<Real,Specs>>& _particles;

    apt::Index<DGrid> _Ib;
    apt::Index<DGrid> _extent;

  public:
    Injector ( const knl::Grid<Real,DGrid>& localgrid,
               const field::Field<Real, 3, DGrid>& Efield,
               const field::Field<Real, 3, DGrid>& Bfield,
               const field::Field<RealJ, 3, DGrid>& Jfield, // J is Jmesh on a replica
               particle::map<particle::array<Real,Specs>>& particles )
      : _grid(localgrid), _Efield(Efield), _Bfield(Bfield), _Jfield(Jfield), _particles(particles) {
      _Ib[0] = pic::ofs::indent[0] - 1;
      _extent[0] = 1;
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
