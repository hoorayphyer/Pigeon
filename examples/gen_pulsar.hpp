#ifndef _GEN_HPP_
#define _GEN_HPP_

#include "apt/numeric.hpp"
#include "apt/index.hpp"

#include "manifold/grid.hpp"

#include "field/field.hpp"

#include "particle/properties.hpp"
#include "particle/map.hpp"
#include "particle/forces.hpp"
#include "particle/scattering.hpp"
#include "particle/particle.hpp"
#include "particle/array.hpp"

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

  inline constexpr const char* project_name = "Pulsar";
  inline constexpr const char* datadir_prefix = "../Data/";

  inline constexpr apt::array<int,DGrid> dims = { 1, 1 };
  inline constexpr apt::array<bool,DGrid> periodic = {false,false};
  inline constexpr int total_timesteps = 1000;
  inline constexpr real_t dt = 0.001;

  constexpr mani::Grid<real_t,DGrid> supergrid
  = {{ { 0.0, std::log(30.0), 128 }, { 0.0, PI, 128 } }};
  inline constexpr int guard = 1;

  inline constexpr int Np = 5;
  inline constexpr real_t epsilon = 1.0 / 5.0;

  constexpr real_t classic_electron_radius () noexcept {
    real_t res = epsilon * epsilon / ( 4 * PI* Np * dt * dt);
    apt::foreach<0,DGrid>( [&res](const auto& g) { res *= g.delta(); }, supergrid );
    return res;
  }
}

namespace pic {
  // TODOL these will become free parameters
  inline constexpr real_t mu0 = 60000.0;
  inline constexpr real_t omega_max = 1.0 / 6.0;

  inline constexpr int spinup_duration = 10.0;

  constexpr real_t omega_spinup ( real_t time ) noexcept {
    return std::min( time / pic::spinup_duration, 1.0 ) * pic::omega_max;
  }
  namespace ofs {
    inline constexpr int magnetic_pole = 2; // 1 for mono-, 2 for di-
    inline constexpr int indent[4] = { 5, 43, guard, guard };
    inline constexpr real_t damping_rate = 0.012;
  }

}

namespace pic {
  inline constexpr int pmpio_num_files = 2;
}

namespace pic :: interval {
  inline constexpr int data_export = 200;
}

// TODOL all the stuff under this {} are meant to be user-specified. Here the pulsar in LogSpherical is used
namespace particle {
  template < typename T, template < typename > class Specs >
  using Pat = typename array< T, Specs >::particle_type; // type of particle from array

  template < typename T, template < typename > class Specs >
  struct CurvatureRadiate : public scat::Channel<T, Specs> {
    using Vec = apt::Vec<T,Specs<T>::Dim>;

    constexpr T calc_Rc ( const Pat<T,Specs>& ptc, const Vec& dp, T dt, const Vec& B ) noexcept {
      // TODOL uniform Rc is used here
      return 1.0;



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

    T (*sample_E_ph)() = nullptr;

    unsigned int abs_charge_x{};
    unsigned int mass_x{};

    T K_thr{};
    T gamma_off{};
    T emission_rate{};

    // return sampled energy if any
    std::optional<T> operator() ( const Pat<T,Specs>& ptc, const Vec& dp, T dt,
                                  const Vec& B, util::Rng<T>& rng ) override {
      T Rc = calc_Rc( ptc, dp, dt, B );
      T gamma = std::sqrt( (mass_x != 0) + apt::sqabs(ptc.p()) );

      if (gamma > gamma_off && gamma > K_thr *  std::cbrt(Rc) && rng.uniform() < emission_rate * dt * gamma  / Rc ) {
        return { std::min( sample_E_ph(), gamma - 1.0 ) };
      } else return {};

    }
  };

  template < typename T, template < typename > class Specs >
  struct MagneticConvert : public scat::Channel<T,Specs> {
    T B_thr{};
    T mfp{};

    std::optional<T> operator() ( const Pat<T,Specs>& photon, const apt::Vec<T,Specs<T>::Dim>& dp, T dt,
                                  const apt::Vec<T,Specs<T>::Dim>& B, util::Rng<T>& rng ) override {
      // prob_mag_conv = dt / mfp_mag_conv
      return ( apt::sqabs(B) > B_thr * B_thr ) && ( rng.uniform() < dt / mfp )  ? std::optional<T>(0.0) : std::nullopt; // NOTE return std::optional(0.0) so as to be treated true in boolean conversion
    }
  };

  template < typename T, template < typename > class Specs >
  struct TwoPhotonCollide : public scat::Channel<T,Specs> {
    // TODO double check this implementation, it is not equivalent because the original one has some sort of gaussian in it. Use Monte Carlo
    // inline T f_x ( T x ) {
    //   // distribution of x*exp(-x^2/2), which peaks at x = 1.
    //   return std::sqrt( -2.0 * std::log(x) );
    // }
    T mfp{};

    std::optional<T> operator() ( const Pat<T,Specs>& photon, const apt::Vec<T,Specs<T>::Dim>& dp, T dt,
                                  const apt::Vec<T,Specs<T>::Dim>& B, util::Rng<T>& rng ) override {
      // prob_mag_conv = dt / mfp_ph_ph
      return ( rng.uniform() < dt / mfp )  ? std::optional<T>(0.0) : std::nullopt;
    }
  };

  // LogSpherical
  namespace force {
    template < typename T, template < typename > class Specs, template < typename, template < typename > class > class Ptc_t >
    void gravity( Ptc_t<T,Specs>& ptc, T dt, const apt::Vec<T,Specs<T>::Dim>& , const apt::Vec<T,Specs<T>::Dim>&, T g  ) noexcept {
      ptc.p()[0] -= g * std::exp( - 2 * ptc.q()[0] ) * dt;
    };

    // when B is strong enough, damp the perpendicular component of momentum
    template < typename T, template < typename > class PtcSpecs, template < typename, template < typename > class > class Ptc_t >
    void landau0( Ptc_t<T,PtcSpecs>& ptc, T dt, const apt::Vec<T,PtcSpecs<T>::Dim>& E, const apt::Vec<T,PtcSpecs<T>::Dim>& B, T B2_thr  ) {
      using Vec = apt::Vec<T,PtcSpecs<T>::Dim>;
      if ( apt::sqabs(B) < B2_thr ) return;

      auto EB2 = apt::dot(E,B);
      EB2 = EB2 * EB2;
      auto B2_E2 = apt::sqabs(B) - apt::sqabs(E);
      // calculate E'^2
      auto Ep2 = 2 * EB2 / ( std::sqrt(B2_E2 * B2_E2 + 4 * EB2) + B2_E2 );
      Vec beta_ExB = apt::cross(E,B) / ( apt::sqabs(B) + Ep2);
      // find B' modulo gamma_ExB
      Vec Bp = B - apt::cross( beta_ExB, E);
      // obtain the momentum with perpendicular components damped
      ptc.p() = Bp * ( apt::dot( ptc.p(), Bp ) / apt::sqabs(Bp) );
      ptc.p() += beta_ExB * std::sqrt( ( 1.0 + apt::sqabs(ptc.p()) ) / ( 1.0 - apt::sqabs(beta_ExB) ) );
    };
  }
}

namespace particle {
  // TODO call this function in particle updater
  template < typename Real >
  void set_up() {
    using namespace pic;
    {
      // TODO initialize properties at a different place
      // TODO move out
      properties[species::electron] = {1,-1,"electron"};
      properties[species::positron] = {1,1,"positron"};
      properties[species::ion] = { 5, 1, "ion"};
      properties[species::photon] = { 0, 0, "photon" };
    }

    {
      using Force = force::Force<real_t,Specs,vParticle>;
      constexpr auto& fgen = force_gen<real_t,Specs,vParticle>;
      constexpr auto* lorentz = force::template lorentz<real_t,Specs,vParticle>;
      Real landau0_B_thr = 0.1 * pic::mu0;
      constexpr auto* landau0 = force::landau0<real_t,Specs,vParticle>;

      constexpr Real gravity_strength = 1.8;
      constexpr auto* gravity = force::gravity<real_t,Specs,vParticle>;
      {
        auto sp = species::electron;
        Force force;
        const auto& prop = properties.at(sp);

        force.add( lorentz, static_cast<Real>(prop.charge_x) / prop.mass_x );
        force.add( gravity, gravity_strength );
        force.add( landau0, landau0_B_thr );

        fgen.Register( sp, force );
      }
      {
        auto sp = species::positron;
        Force force;
        const auto& prop = properties.at(sp);

        force.add( lorentz, static_cast<Real>(prop.charge_x) / prop.mass_x );
        force.add( gravity, gravity_strength );
        force.add( landau0, landau0_B_thr );

        fgen.Register( sp, force );
      }
      {
        auto sp = species::ion;
        Force force;
        const auto& prop = properties.at(sp);

        force.add( lorentz, static_cast<Real>(prop.charge_x) / prop.mass_x );
        force.add( gravity, gravity_strength );
        // force.add( landau0, landau0_B_thr );

        fgen.Register( sp, force );
      }
    }


    {
      scat::RadiationFromCharges<false,real_t,Specs> ep_scat;

      {
        struct InZone : public scat::Eligible<real_t,Specs> {
          bool operator() ( const Pat<real_t,Specs>& ptc ) override {
            return ptc.q()[0] < std::log(9.0);
          }
        };
        ep_scat.add(InZone{});
      }

      {
        CurvatureRadiate<real_t,Specs> cr;
        cr.abs_charge_x = 1u;
        cr.mass_x = 1u;
        cr.K_thr = 80.0;
        cr.gamma_off = 15.0;
        cr.emission_rate = 0.25;
        cr.sample_E_ph = []() noexcept { return 14.0; };

        ep_scat.add(cr);
      }

      scat_gen<real_t,Specs>.Register( species::electron, ep_scat );
      scat_gen<real_t,Specs>.Register( species::positron, ep_scat );
    }

    {
      scat::PhotonPairProduction<real_t,Specs> photon_scat;
      // Photons are free to roam across all domain. They may produce pairs outside light cylinder
      {
        MagneticConvert<real_t,Specs> mc;
        mc.B_thr = pic::mu0 / 10.0 ;
        mc.mfp = 0.2;
        photon_scat.add(mc);
      }
      {
        TwoPhotonCollide<real_t,Specs> tc;
        tc.mfp = 5.0;
        photon_scat.add(tc);
      }

      scat_gen<real_t,Specs>.Register( species::photon, photon_scat );
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

    Real B_r_over_mu0 ( Real logr, Real theta ) noexcept {
      return 2.0 * std::cos(theta) * std::exp(-3.0 * logr);
    };

    Real B_th_over_mu0 ( Real logr, Real theta ) noexcept {
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
        _Bfield[0](I) = _mu0 * B_r_over_mu0( _grid[0].absc(I[0], _Bfield[0].offset()[0]), _grid[1].absc(I[1], _Bfield[0].offset()[1]) );
        _Bfield[1](I) = _mu0 * B_th_over_mu0( _grid[0].absc(I[0], _Bfield[1].offset()[0]), _grid[1].absc(I[1], _Bfield[1].offset()[1]) );
      }
    }

    constexpr int initial_timestep() noexcept { return 0; }
  };

  template < int DGrid,
             typename Real,
             template < typename > class Specs,
             typename RealJ >
  struct FieldBC_Rotating_Conductor {
  private:
    const mani::Grid<Real,DGrid>& _grid;
    field::Field<Real, 3, DGrid>& _Efield;
    field::Field<Real, 3, DGrid>& _Bfield;

    apt::Index<DGrid> _Ib;
    apt::Index<DGrid> _extent;

    const Real _mu0 = pic::mu0;

    Real B_r_over_mu0 ( Real logr, Real theta ) noexcept {
      return 2.0 * std::cos(theta) * std::exp(-3.0 * logr);
    };

    Real B_th_over_mu0 ( Real logr, Real theta ) noexcept {
      return std::sin(theta) * std::exp(-3.0 * logr);
    };

    // find out E by E = -( Omega x r ) x B
    Real E_r_over_mu0omega ( Real logr, Real theta ) noexcept {
      auto sin_t = std::sin(theta);
      return std::exp( -2.0 * logr ) * sin_t * sin_t;
    };

    Real E_th_over_mu0omega ( Real logr, Real theta ) noexcept {
      return - std::exp( -2.0 * logr ) * std::sin( 2.0 * theta );
    };

  public:
    FieldBC_Rotating_Conductor ( const mani::Grid<Real,DGrid>& localgrid,
                                 field::Field<Real, 3, DGrid>& Efield,
                                 field::Field<Real, 3, DGrid>& Bfield,
                                 const field::Field<RealJ, 3, DGrid>& Jfield, // J is Jmesh on a replica
                                 const particle::map<particle::array<Real,Specs>>& particles )
      : _grid(localgrid), _Efield(Efield), _Bfield(Bfield) {
          _Ib[0] = 0;
          _extent[0] = pic::ofs::indent[0];
          apt::tie(_Ib[1], _extent[1]) = gtl( {0.0, PI}, localgrid[1] );
    }

    void operator() ( int timestep, Real dt ) {
      const auto omega = pic::omega_spinup( timestep * dt );
      for ( auto I : apt::Block(_extent) ) {
        I += _Ib;
        // TODO deal with discontinuous and continuous variables
        // TODO piecing together different patches
        _Bfield[0](I) = _mu0 * B_r_over_mu0( _grid[0].absc(I[0], _Bfield[0].offset()[0]), _grid[1].absc(I[1], _Bfield[0].offset()[1]) );
        _Bfield[1](I) = _mu0 * B_th_over_mu0( _grid[0].absc(I[0], _Bfield[1].offset()[0]), _grid[1].absc(I[1], _Bfield[1].offset()[1]) );
        _Bfield[2](I) = 0.0;

        _Efield[0](I) = _mu0 * omega * E_r_over_mu0omega( _grid[0].absc(I[0], _Efield[0].offset()[0]), _grid[1].absc(I[1], _Efield[0].offset()[1]) );
        _Efield[1](I) = _mu0 * omega * E_th_over_mu0omega( _grid[0].absc(I[0], _Efield[1].offset()[0]), _grid[1].absc(I[1], _Efield[1].offset()[1]) );
        _Efield[2](I) = 0.0;
      }
    }
  };

  template < int DGrid,
             typename Real,
             template < typename > class Specs,
             typename RealJ >
  struct FieldBC_Axis {
  private:
    const mani::Grid<Real,DGrid>& _grid;
    field::Field<Real, 3, DGrid>& _Efield;
    field::Field<Real, 3, DGrid>& _Bfield;

    const int axis_dir = 1;
    bool _is_at_axis_lower = false;
    bool _is_at_axis_upper = false;

  public:
    FieldBC_Axis ( const mani::Grid<Real,DGrid>& localgrid,
                        field::Field<Real, 3, DGrid>& Efield,
                        field::Field<Real, 3, DGrid>& Bfield,
                        const field::Field<RealJ, 3, DGrid>& Jfield,
                        const particle::map<particle::array<Real,Specs>>& particles )
      : _grid(localgrid), _Efield(Efield), _Bfield(Bfield) {
      _is_at_axis_lower = std::abs( localgrid[axis_dir].lower() - 0.0 ) < localgrid[axis_dir].delta();
      _is_at_axis_upper = std::abs( localgrid[axis_dir].upper() - PI ) < localgrid[axis_dir].delta();
    }

    void operator() () {
      // TODO Guard cells values are needed when doing interpolating E and B
      // E_theta, B_r, B_phi are on the axis. All but B_r should be set to zero
      const auto& mesh = _Efield.mesh();
      if ( _is_at_axis_lower ) {
        int n = 0;
        for ( const auto& trI : mesh.project(1, {}, mesh.extent() ) ) {
          _Efield[1][trI | n] = 0.0;
          _Bfield[2][trI | n] = 0.0;
        }
      }

      if ( _is_at_axis_upper ) {
        int n = _grid[1].dim();
        for ( const auto& trI : mesh.project(1, {}, mesh.extent() ) ) {
          _Efield[1][trI | n] = 0.0;
          _Bfield[2][trI | n] = 0.0;
        }
      }
    }
  };

  template < int DGrid,
             typename Real,
             template < typename > class Specs,
             typename RealJ >
  struct FieldBC_FoldBackJ {
  private:
    static_assert(DGrid==2);
    const mani::Grid<Real,DGrid>& _grid;
    field::Field<RealJ, 3, DGrid>& _Jmesh;

    const int axis_dir = 1;
    bool _is_at_axis_lower = false;
    bool _is_at_axis_upper = false;

  public:
    FieldBC_FoldBackJ ( const mani::Grid<Real,DGrid>& localgrid,
                        const field::Field<Real, 3, DGrid>& Efield,
                        const field::Field<Real, 3, DGrid>& Bfield,
                        field::Field<RealJ, 3, DGrid>& Jfield,
                        const particle::map<particle::array<Real,Specs>>& particles )
      : _grid(localgrid), _Jmesh(Jfield) {
      _is_at_axis_lower = std::abs( localgrid[axis_dir].lower() - 0.0 ) < localgrid[axis_dir].delta();
      _is_at_axis_upper = std::abs( localgrid[axis_dir].upper() - PI ) < localgrid[axis_dir].delta();
    }

    void operator() () {
      const auto& mesh = _Jmesh.mesh();
      const int guard = mesh.guard();

      if ( _is_at_axis_lower ) {
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
      }

      if ( _is_at_axis_upper ) {
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
      constexpr Real j_reg_x = 2.0;
      constexpr Real Ninj = 5.0;

      constexpr auto posion = species::ion;
      constexpr auto negaon = species::electron;

      Real omega = pic::omega_spinup( timestep * dt );

      int charge_x = particle::properties.at(posion).charge_x;

      // the idea is that (timestep + 1) * N_inj * profile represents the accummulated number of injected pairs through the specified timestep. NOTE the timestep is shifted by one to reflect the actual times the inject is called, including this time.
      // NOTE _Jfield and Jmesh also differs in their underlying mesh guard cells
      auto profile_inj =
        [&Ninj,&timestep] ( Real theta ) noexcept {
          Real inj_num_base = Ninj * 0.5 * std::abs( std::sin( 2 * theta ) );
          return static_cast<int>( (timestep + 1) * inj_num_base ) - static_cast<int>( timestep * inj_num_base );
        };

      auto j_reg_inj =
        [&grid=_grid, charge_x, &j_reg_x]( Real J ) noexcept {
          static auto factor = j_reg_x * grid[0].delta() * grid[1].delta() / charge_x;
          return J * factor;
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

        int num = std::max<Real>( profile_inj(q[1]), j_reg_inj(apt::abs(J)) );

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
