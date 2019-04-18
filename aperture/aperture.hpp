#ifndef  _APERTURE_HPP_
#define  _APERTURE_HPP_

#include "abstract_field_updater.hpp"
#include "kernel/grid.hpp"

#include "old_field_solver/old_field_solver_adapter.hpp"
#include "particle_updater.hpp"
#include "field/communication.hpp"

#include "dye/dynamic_balance.hpp"

#include "particle/migration.hpp"

#include <memory>

#include "apt/vec.hpp"
#include "apt/numeric.hpp"
#include "particle/particle.hpp"

namespace aperture {
  constexpr long double PI = std::acos(-1.0l);

  template < typename Real, int DGrid >
  apt::pair< apt::Index<DGrid> > gtl ( const apt::array<apt::pair<Real>,DGrid>& range,
                                       const knl::Grid<Real,DGrid>& localgrid ) noexcept {
    apt::Index<DGrid> Ib;
    apt::Index<DGrid> extent;
    for ( int i = 0; i < DGrid; ++i ) {
      Ib[i] = std::max<int>( 0, ( range[i][LFT] - localgrid[i].lower() ) / localgrid[i].delta() );
      extent[i] = std::min<int>( range[i][RGT], localgrid[i].dim() ) - Ib[i];
    }
    return { Ib, extent };
  }

  template < typename Real, int DGrid, typename state_t, typename RealJ >
  struct InitialConditionDipole {
  private:
    const knl::Grid<Real,DGrid>& _grid;
    field::Field<Real, 3, DGrid>& _Bfield;

    apt::Index<DGrid> _Ib;
    apt::Index<DGrid> _extent;

    const Real _mu0 = 60000.0;

    Real B_r_over_mu0 ( Real logr, Real theta ) noexcept {
      return 2.0 * std::cos(theta) * std::exp(-3.0 * logr);
    };

    Real B_th_over_mu0 ( Real logr, Real theta ) noexcept {
      return std::sin(theta) * std::exp(-3.0 * logr);
    };

  public:
    InitialConditionDipole ( const knl::Grid<Real,DGrid>& localgrid,
                             const field::Field<Real, 3, DGrid>& Efield,
                             field::Field<Real, 3, DGrid>& Bfield,
                             const field::Field<RealJ, 3, DGrid>& Jfield, // J is Jmesh on a replica
                             const particle::map<particle::array<Real,3,state_t>>& particles )
      : _grid(localgrid), _Bfield(Bfield) {
      constexpr apt::array<apt::pair<Real>,DGrid> global_range = {{ {1.0, std::log(30.0)}, {0, PI} }};
      apt::tie(_Ib, _extent) = gtl( global_range, localgrid );
    }

    void operator() () {
      for ( auto I : apt::Block(_extent) ) {
        I += _Ib;
        _Bfield[0](I) = _mu0 * B_r_over_mu0( _grid[0].absc(I[0], _Bfield[0].offset()[0]), _grid[1].absc(I[1], _Bfield[0].offset()[1]) );
        _Bfield[1](I) = _mu0 * B_th_over_mu0( _grid[0].absc(I[0], _Bfield[1].offset()[0]), _grid[1].absc(I[1], _Bfield[1].offset()[1]) );
      }
    }
  };

  template < typename Real, int DGrid, typename state_t, typename RealJ >
  struct FieldBC_Rotating_Conductor {
  private:
    const knl::Grid<Real,DGrid>& _grid;
    field::Field<Real, 3, DGrid>& _Efield;
    field::Field<Real, 3, DGrid>& _Bfield;

    apt::Index<DGrid> _Ib;
    apt::Index<DGrid> _extent;

    const Real _mu0 = 60000.0;

    Real omega_spinup ( Real time ) noexcept {
      return std::min( time / 10.0, 1.0 ) / 6.0;
    }

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
    FieldBC_Rotating_Conductor ( const knl::Grid<Real,DGrid>& localgrid,
                                 field::Field<Real, 3, DGrid>& Efield,
                                 field::Field<Real, 3, DGrid>& Bfield,
                                 const field::Field<RealJ, 3, DGrid>& Jfield, // J is Jmesh on a replica
                                 const particle::map<particle::array<Real,3,state_t>>& particles )
      : _grid(localgrid), _Efield(Efield), _Bfield(Bfield) {
      constexpr apt::array<apt::pair<Real>,DGrid> global_range = {{ {0.0, std::log(1.01)}, {0, PI} }};
      apt::tie(_Ib, _extent) = gtl( global_range, localgrid );
    }

    void operator() ( int timestep, Real dt ) {
      const auto omega = omega_spinup( timestep * dt );
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

  template < bool IsLower, typename Real, int DGrid, typename state_t, typename RealJ >
  struct FieldBC_Axis {
  private:
    const knl::Grid<Real,DGrid>& _grid;
    field::Field<Real, 3, DGrid>& _Efield;
    field::Field<Real, 3, DGrid>& _Bfield;

    bool _is_at_axis = false;

  public:
    FieldBC_Axis ( const knl::Grid<Real,DGrid>& localgrid,
                        field::Field<Real, 3, DGrid>& Efield,
                        field::Field<Real, 3, DGrid>& Bfield,
                        const field::Field<RealJ, 3, DGrid>& Jfield,
                        const particle::map<particle::array<Real,3,state_t>>& particles )
      : _grid(localgrid), _Efield(Efield), _Bfield(Bfield) {
      if constexpr ( IsLower )
                     _is_at_axis = std::abs( localgrid[1].lower() - 0.0 ) < localgrid[1].delta();
      else
        _is_at_axis = std::abs( localgrid[1].upper() - PI ) < localgrid[1].delta();
    }

    void operator() () {
      // TODO We don't need to do anything to the guard cells right?
      if ( !_is_at_axis ) return;
      // E_theta, B_r, B_phi are on the axis. All but B_r should be set to zero
      const auto& mesh = _Efield.mesh();
      int n = IsLower ? 0 : _grid[1].dim();
      for ( const auto& trI : mesh.project(1, {}, mesh.extent() ) ) {
        _Efield[1][trI | n] = 0.0;
        _Bfield[2][trI | n] = 0.0;
      }
    }
  };

  template < bool IsLower, typename Real, int DGrid, typename state_t, typename RealJ >
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
                        const particle::map<particle::array<Real,3,state_t>>& particles )
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
  template < typename Real, int DPtc, int DGrid, typename state_t, typename RealJ >
  struct Injector {
  private:
    const knl::Grid<Real,DGrid>& _grid;
    const field::Field<Real, 3, DGrid>& _Efield;
    const field::Field<Real, 3, DGrid>& _Bfield;
    const field::Field<RealJ, 3, DGrid>& _Jfield; // J is Jmesh on a replica
    particle::map<particle::array<Real,3,state_t>>& _particles;

    apt::Index<DGrid> _Ib;
    apt::Index<DGrid> _extent;

    Real omega_spinup ( Real time ) noexcept {
      return std::min( time / 10.0, 1.0 ) / 6.0;
    }

  public:
    Injector ( const knl::Grid<Real,DGrid>& localgrid,
               const field::Field<Real, 3, DGrid>& Efield,
               const field::Field<Real, 3, DGrid>& Bfield,
               const field::Field<RealJ, 3, DGrid>& Jfield, // J is Jmesh on a replica
               particle::map<particle::array<Real,3,state_t>>& particles )
      : _grid(localgrid), _Efield(Efield), _Bfield(Bfield), _Jfield(Jfield), _particles(particles) {
      constexpr apt::array<apt::pair<Real>,DGrid> global_range = {{ {std::log(1.01), std::log(1.02)}, {0, PI} }};
      apt::tie(_Ib, _extent) = gtl( global_range, localgrid );
    }

    void operator() ( int timestep, Real dt, util::Rng<Real>& rng ) {
      using namespace particle;

      constexpr Real v_th = 0.3;
      constexpr Real j_reg_x = 2.0;
      constexpr int Ninj = 10;

      constexpr auto posion = species::positron;
      constexpr auto negaon = species::electron;

      Real omega = omega_spinup( timestep * dt );

      int charge_x = 1.0;

      // the idea is that (timestep + 1) * N_inj * profile represents the accummulated number of injected pairs through the specified timestep. NOTE the timestep is shifted by one to reflect the actual times the inject is called, including this time.
      // NOTE _Jfield and Jmesh also differs in their underlying mesh guard cells
      auto profile_inj =
        [&Ninj,&timestep] ( Real theta ) noexcept {
          Real inj_num_base = Ninj * 0.5 * std::abs( std::sin( 2 * theta ) );
          return static_cast<int>( (timestep + 1) * inj_num_base ) - static_cast<int>( timestep * inj_num_base );
        };

      auto j_reg_inj =
        [&grid=_grid, charge_x, &j_reg_x]( Real J ) noexcept {
          static auto factor = j_reg_x * grid[0].delta() * grid[1].delta() / std::labs(charge_x);
          return J * factor;
        };

      auto itr_po = std::back_inserter(_particles[posion]);
      auto itr_ne = std::back_inserter(_particles[negaon]);

      for ( auto I : apt::Block(_extent) ) {
        I += _Ib;
        apt::Vec<Real, DPtc> q{};
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

        apt::Vec<Real, DPtc> p{};
        p[2] = omega * std::exp(q[0]) * std::sin(q[1]); // corotating

        for ( int n = 0; n < num; ++n ) {
          auto q_ptc = q;
          for ( int i = 0; i < DGrid; ++i )
            q_ptc[i] += rng.uniform();
          auto p_ptc = p;
          p_ptc += nB * rng.gaussian( 0.0, v_th );
          *(itr_ne++) = Particle<Real,DPtc,state_t>( q_ptc, p_ptc, negaon );
          *(itr_po++) = Particle<Real,DPtc,state_t>( std::move(q_ptc), std::move(p_ptc), posion );
        }
      }
    }
  };
}

// TODOL merge this into Aperture
#include "traits.hpp"
namespace aperture {
  using namespace traits;
  template < typename Real, int DGrid, typename state_t >
  using ParticleUpdater = particle::ParticleUpdater< Real, DGrid, 3, state_t, ShapeF, real_j_t, coordinate_system >;
}

namespace aperture {
  template < typename Real, int DGrid, typename state_t, typename RealJ >
  struct Aperture {
  private:
    const knl::Grid< Real, DGrid >& _supergrid;
    const int _guard;
    std::optional<mpi::CartComm> _cart_opt;
    util::Rng<Real> _rng;

    knl::Grid< Real, DGrid > _grid;
    std::optional<dye::Ensemble<DGrid>> _ens_opt;
    apt::array< apt::pair<Real>, DGrid > _borders;

    field::Field<Real, 3, DGrid> _E;
    field::Field<Real, 3, DGrid> _B;
    field::Field<RealJ, 3, DGrid> _J;
    particle::map<particle::array<Real,3,state_t>> _particles;

    std::unique_ptr<AbstractFieldUpdater<Real,DGrid,RealJ>> _field_update;
    std::unique_ptr<ParticleUpdater<Real, DGrid, state_t>> _ptc_update;

    std::vector<particle::cParticle<Real, DPtc, state_t>> _migrators;

    FieldBC_FoldBackJ<true, Real, DGrid, state_t, RealJ> _fbj_lower;
    FieldBC_FoldBackJ<false, Real, DGrid, state_t, RealJ> _fbj_upper;
    Injector<Real, 3, DGrid, state_t, RealJ> _injector;

    Real _unit_e;

    // ScalarField<Scalar> pairCreationEvents; // record the number of pair creation events in each cell.
    // PairCreationTracker pairCreationTracker;

    constexpr auto all_but( field::offset_t all, int but_comp ) noexcept {
      apt::array<field::offset_t, DGrid> res;
      apt::foreach<0,DGrid>( [&all]( auto& ofs ) { ofs = all; }, res );
      if ( but_comp < DGrid ) res[but_comp] = !all;
      return res;
    }


    void refresh( const dye::Ensemble<DGrid>& ens, double unit_e ) {

      apt::Index<DGrid> bulk_dims;
      for ( int i = 0; i < DGrid; ++i ) {
        int dim = _supergrid[i].dim() / ens.cart_dims[i];
        _grid[i] = _supergrid[i];
        _grid[i].clip( ens.cart_coords[i] * dim, dim );
        // TODOL cart_dim = 1 and periodic
        _borders[i] = { _grid[i].lower(), _grid[i].upper() };
        bulk_dims[i] = _grid[i].dim();
      }

      {
        _E = {{ bulk_dims, _guard }};
        _B = {{ bulk_dims, _guard }};
        // NOTE minimum number of guards of J on one side is ( supp + 1 ) / 2 + 1
        _J = {{ bulk_dims, ( ShapeF::support() + 3 ) / 2 }};

        for( int i = 0; i < 3; ++i ) {
          _E.set_offset( i, all_but( MIDWAY, i ) );
          _B.set_offset( i, all_but( INSITU, i ) );
          _J.set_offset( i, all_but( MIDWAY, i ) );
        }
        _E.reset();
        _B.reset();
        _J.reset();
      }

      if ( _cart_opt )
        _field_update.reset(new ofs::OldFieldUpdater<>( unit_e, *_cart_opt, _grid, ens.is_at_boundary(), _guard ) );
      _ptc_update.reset(new ParticleUpdater<Real, DGrid, state_t>( _grid, _rng ) );
    }

  public:
    Aperture( const knl::Grid< Real, DGrid >& supergrid, const std::optional<mpi::CartComm>& cart_opt, int guard, util::Rng<Real> rng, double unit_e )
      : _supergrid(supergrid), _guard(guard), _cart_opt(cart_opt), _rng(std::move(rng)),
        _injector{ _grid, _E, _B, _J, _particles },
        _fbj_lower{ _grid, _E, _B, _J, _particles },
        _fbj_upper{ _grid, _E, _B, _J, _particles },
        _unit_e(unit_e) {
      _grid = supergrid;
      _ens_opt = dye::create_ensemble<DGrid>(cart_opt);
      if ( !_ens_opt ) return;

      const auto& ens = *_ens_opt;
      refresh(ens, unit_e);

      InitialConditionDipole ic( _grid, _E, _B, _J, _particles );
      ic();
    }

    void evolve( int timestep, Real dt ) {
      if ( _ens_opt ) {
        const auto& ens = *_ens_opt;

        // TODOL Opimize communication. Use persistent and buffer?
        for ( int i = 0; i < 3; ++i ) {
          auto& buffer = _J[i].data();
          ens.intra.template reduce<mpi::IN_PLACE>( mpi::by::SUM, ens.chief, buffer.data(), buffer.size() );
        }

        if ( _cart_opt ) {
          field::merge_guard_cells_into_bulk( _J, *_cart_opt );
          (*_field_update)(_E, _B, _J, dt, timestep);
        }
        // TODOL reduce number of communications?
        for ( int i = 0; i < 3; ++i )
          ens.intra.broadcast( ens.chief, _E[i].data().data(), _E[i].data().size() );
        for ( int i = 0; i < 3; ++i )
          ens.intra.broadcast( ens.chief, _B[i].data().data(), _B[i].data().size() );

        _J.reset();
        // if ( false )
        //   sort_particles();
        (*_ptc_update) ( _particles, _J, _E, _B, dt, timestep );

        _fbj_lower();
        _fbj_upper();

        { // migration
          auto migrate_dir =
            []( auto q, auto lb, auto ub ) noexcept {
              return ( q >= lb ) + ( q > ub );
            };

          for ( auto&[ sp, ptcs ] : _particles ) {
            for ( auto ptc : ptcs ) { // TODOL semantics
              char mig_dir = 0;
              for ( int i = 0; i < DGrid; ++i ) {
                mig_dir += migrate_dir( ptc.q()[i], _borders[i][LFT], _borders[i][RGT] ) * apt::pow3(i);
              }

              if ( mig_dir != ( apt::pow3(DGrid) - 1 ) / 2 ) {
                _migrators.emplace_back(std::move(ptc));
                _migrators.back().extra() = mig_dir;
              }
            }
          }

          particle::migrate( _migrators, ens.inter, timestep );
          for ( auto&& ptc : _migrators ) {
            _particles[ptc.template get<particle::species>()].push_back( std::move(ptc) );
          }
          _migrators.resize(0);
        }

        _injector( timestep, dt, _rng );
      }


      // TODO NOTE injection and annihilation are more like free sources and sinks of particles users control, in other words they fit the notion of particle boundary conditions. Do them outside of Updater. In contrast, pair creation is a physical process we model, so it is done inside the updater.
      // inject_particles(); // TODO only coordinate space current is needed in implementing current regulated injection // TODO compatible with ensemble??
      // TODOL annihilation will affect deposition // NOTE one can deposit in the end
      // annihilate_mark_pairs( );

      if (false) {
        // TODO check idle?
        // TODO export_data();
        if ( false ) {
          // TODO has a few hyper parameters
          // TODO touch create is not multinode safe even buffer is used
          std::optional<int> old_label;
          if ( _ens_opt ) old_label.emplace(_ens_opt->label());

          dynamic_load_balance( _particles, _ens_opt, _cart_opt, 100000 );

          std::optional<int> new_label;
          if ( _ens_opt ) new_label.emplace(_ens_opt->label());
          if ( old_label != new_label ) refresh(*_ens_opt, _unit_e);
        }
      }

      // TODOL
      // if (false)
      //   save_snapshot();

      // TODOL
      // if (false)
      //   save_tracing();
    }
  };
}

#endif
