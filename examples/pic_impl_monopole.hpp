#ifndef _PIC_IMPL_HPP_
#define _PIC_IMPL_HPP_

#include "../examples/pic_prior_impl.hpp"
#include "apt/numeric.hpp"
#include "metric/log_spherical.hpp"
#include "field/log_spherical_solver/updater.hpp"

#include "pic/forces/gravity.hpp"
#include "io/exportee_by_function.hpp"
#include "io/data_exporter.hpp"
#include "io/exportee.hpp"
#include "msh/mesh_shape_interplay.hpp"
#include <cassert>


namespace pic {
  using Metric = metric::LogSpherical<real_t>;

  inline constexpr apt::array<int,DGrid> dims = { 1, 1 };
  inline constexpr apt::array<bool,DGrid> periodic = {false,false};

  constexpr Grid supergrid
  = {{ { 0.0, std::log(30.0), 128 }, { 0.0, 180.0_deg, 128 } }};

  real_t Omega;
  real_t dt;
  real_t mu;
  real_t wpic2;

  real_t gravity_strength;
  real_t atm_x;
  real_t v_th;

  real_t damping_thickness;
  real_t damping_rate;
  real_t spinup_time;

  std::string project_name;
  std::string datadir_prefix;
  int total_timesteps;

  void load_configuration(const ConfFile_t& conf) {
    safe_set(dt, conf["dt"]);
    safe_set(Omega, conf["Omega"]);
    safe_set(mu, conf["mu"]);
    safe_set(wpic2, conf["Np"]); wpic2 = 2 * Omega * mu / wpic2;

    safe_set(gravity_strength, conf["forces"]["gravity"]);

    safe_set(v_th, conf["atmosphere"]["v_th"]);
    safe_set(atm_x, conf["atmosphere"]["multiplicity"]);

    safe_set(damping_thickness, conf["damping"]["thickness"]);
    safe_set(damping_rate, conf["damping"]["rate"]);
    safe_set(spinup_time, conf["spinup_time"]);

    project_name = conf["project_name"].value_or("Unnamed"sv);
    datadir_prefix = conf["datadir_prefix"].value_or("../Data/"sv);
    total_timesteps = conf["total_timesteps"].value_or(100);

    print_timestep_to_stdout_interval = conf["plans"]["print_timestep_to_stdout_interval"].value_or(100);

    sort_ptcs_plan.on = conf["plans"]["sort"]["on"].value_or(true);
    sort_ptcs_plan.start = conf["plans"]["sort"]["start"].value_or(0);
    safe_set(sort_ptcs_plan.interval, conf["plans"]["sort"]["interval"]);

    export_plan.on = conf["plans"]["export"]["on"].value_or(true);
    export_plan.start = conf["plans"]["export"]["start"].value_or(0);
    safe_set(export_plan.interval, conf["plans"]["export"]["interval"]);
    export_plan.num_files = conf["plans"]["export"]["num_files"].value_or(1);
    export_plan.downsample_ratio = conf["plans"]["export"]["downsample_ratio"].value_or(1);

    checkpoint_plan.on = conf["plans"]["checkpoint"]["on"].value_or(false);
    checkpoint_plan.start = conf["plans"]["checkpoint"]["start"].value_or(1);
    safe_set(checkpoint_plan.interval, conf["plans"]["checkpoint"]["interval"]);
    checkpoint_plan.num_files = conf["plans"]["checkpoint"]["num_files"].value_or(1);
    checkpoint_plan.max_num_checkpoints = conf["plans"]["checkpoint"]["max_num_checkpoints"].value_or(1);

    load_balance_plan.on = conf["plans"]["load_balance"]["on"].value_or(true);
    load_balance_plan.start = conf["plans"]["load_balance"]["start"].value_or(0);
    safe_set(load_balance_plan.interval, conf["plans"]["load_balance"]["interval"]);
    load_balance_plan.target_load = conf["plans"]["load_balance"]["target_load"].value_or(100000);

    profiling_plan.on = conf["plans"]["profiling"]["on"].value_or(false);
    profiling_plan.start = conf["plans"]["profiling"]["start"].value_or(0);
    safe_set(profiling_plan.interval, conf["plans"]["profiling"]["interval"]);
    if ( auto n = conf["plans"]["profiling"]["max_entries"].value<int64_t>() ) {
      profiling_plan.max_entries = {*n};
    }

    vitals_plan.on = conf["plans"]["vitals"]["on"].value_or(true);
    vitals_plan.start = conf["plans"]["vitals"]["start"].value_or(0);
    safe_set(vitals_plan.interval, conf["plans"]["vitals"]["interval"]);
  }

  real_t r_e() {
    real_t res = wpic2 / ( 4.0 * std::acos(-1.0l));
    res *= apt::dV(supergrid);
    return res;
  }
}

namespace pic {
  constexpr int star_interior = 5;

  constexpr int field_op_inv_precision = 4;
  constexpr int myguard = std::max(
                                   ::field::LogSphericalSolver<real_t, DGrid, real_j_t>::min_guard(field_op_inv_precision),
                                   (pic::ShapeF::support() + 3) / 2); // NOTE minimum number of guards of J on one side is ( supp + 3 ) / 2

  real_t omega_spinup ( real_t time ) noexcept {
    return std::min<real_t>( time / spinup_time, 1.0_r ) * Omega;
  }

  real_t B_r_star( real_t lnr, real_t, real_t , real_t time ) noexcept {
    return pic::mu * std::exp(-2.0_r * lnr);
  }
  real_t B_theta_star( real_t lnr, real_t, real_t , real_t time ) noexcept {
    return 0;
  }
  constexpr real_t B_phi_star( real_t lnr, real_t, real_t , real_t time ) noexcept {
    return 0;
  }

  real_t E_r_star( real_t lnr, real_t theta, real_t , real_t time ) noexcept {
    return 0;
  }
  real_t E_theta_star( real_t lnr, real_t theta, real_t , real_t time ) noexcept {
    return - pic::mu * omega_spinup(time) * std::exp(- lnr ) * std::sin( theta );
  }
  constexpr real_t E_phi_star( real_t lnr, real_t theta, real_t , real_t time ) noexcept {
    return 0;
  }

  template < typename T, typename F >
  void axissymmetrize( ::field::Component<T,DGrid,false> comp, // TODOL semantics on comp
                       const F& f, // has interface: void (*f)( real_t& val_guard, real_t& val_bulk ),
                       Index Ib, Index Ie, bool is_upper ) {
    static_assert(DGrid==2);
    constexpr int AxisDir = 1;

    int mirror_sum = (comp.offset()[AxisDir] == MIDWAY ) ? -1 : 0;
    mirror_sum += is_upper ? 2 * Ib[AxisDir] : 2 * (Ie[AxisDir] - 1);

    for ( const auto& trI : apt::project_out(AxisDir,Ib,Ie) ) {
      for ( apt::Longidx n (AxisDir, Ib[AxisDir]); n < Ie[AxisDir]; ++n ) {
        f( comp(trI + n), comp(trI + (mirror_sum - n)) );
      }
    }
  }
}

namespace pic {
  struct RTD { // runtime data
  public:
    static RTD& data() {
      static RTD r;
      return r;
    }

    map<real_t> N_scat {};

    std::optional<map<JField>> Jsp; // current by species

    void init( const map<Properties>& properties, const Grid& localgrid ) {
      for ( auto sp : properties )
        N_scat.insert( sp, 0 );

      if (false) {// activate Jsp
        Jsp.emplace(std::remove_reference_t<decltype(*Jsp)>{});
        for ( auto sp : properties )
          (*Jsp).insert(sp, {});
      }
    };

  private:
    RTD() = default;
    RTD(const RTD&);
    RTD(RTD&&) noexcept;
    RTD& operator=( const RTD& );
    RTD& operator=( RTD&& ) noexcept;
    ~RTD() = default;
  };
}

namespace pic {
  auto set_up_field_actions() {
    std::vector<std::unique_ptr<FieldAction>> fus;
    namespace range = apt::range;

    ::field::LogSphericalSolver<real_t,DGrid,real_j_t> fu;
    {
      const int guard = myguard;
      fu.setName("LogSphericalSolver");
      fu[0] = { star_interior, pic::supergrid[0].dim(), guard };
      fu[1] = { 0, pic::supergrid[1].dim()+1, guard };

      fu.set_fourpi(4.0 * std::acos(-1.0l) * pic::r_e());
      fu.set_alpha(1.0);
      fu.set_op_inv_precision(field_op_inv_precision);
      fu.set_surface(supergrid[0].absc(star_interior));
      fu.set_outer(supergrid[0].upper());
    }

    // NOTE only implemented for UPPER boundary
    struct RotatingConductor : public FieldAction {
    private:
      apt::array<real_t(*)(real_t,real_t,real_t,real_t t),3> _E_cond;
      apt::array<real_t(*)(real_t,real_t,real_t,real_t t),3> _B_cond;

    public:
      auto& set_E_cond( apt::array<real_t(*)(real_t,real_t,real_t,real_t t),3> E ) {_E_cond = E; return *this;}
      auto& set_B_cond( apt::array<real_t(*)(real_t,real_t,real_t,real_t t),3> B ) {_B_cond = B; return *this;}

      virtual RotatingConductor* Clone() const { return new RotatingConductor(*this); }

      virtual void operator() ( Field<3>& E, Field<3>& B, JField&,
                                const Grid& grid, const mpi::CartComm&,
                                int timestep, real_t dt) const override {
        auto impose =
          [this, &grid, time=timestep*dt]( auto& F, const auto& F_cond, int comp ) {
            const auto& ofs = F[comp].offset();
            const auto& m = F.mesh();

            for ( int j = range::begin(*this,1); j < range::end(*this,1); ++j ) {
              int li_j = m.linear_index(1, j);
              real_t th = grid[1].absc(j, ofs[1]*0.5);
              for ( int i = range::begin(*this,0); i < range::end(*this,0); ++i ) {
                real_t lnr = grid[0].absc(i, ofs[0]*0.5);
                int li = li_j + i * m.stride(0);
                F[comp][li] = F_cond[comp](lnr,th,0,time);
              }
            }
          };
        for ( int C = 0; C < 3; ++C )
          impose(E,_E_cond,C);
        for ( int C = 0; C < 3; ++C )
          impose(B,_B_cond,C);
      }
    } fu_cond;
    {
      auto &fu = fu_cond;

      fu.setName("RotatingConductor");
      fu[0] = {-myguard, star_interior};
      fu[1] = {0, pic::supergrid[1].dim() + 1};

      fu.set_E_cond({E_r_star, E_theta_star, E_phi_star});
      fu.set_B_cond({B_r_star, B_theta_star, B_phi_star});
    }

    // NOTE only implemented for UPPER boundary
    struct DampingLayer : public FieldAction {
    private:
      apt::array<real_t(*)(real_t,real_t,real_t,real_t t),3> _E_bg;
      apt::array<real_t(*)(real_t,real_t,real_t,real_t t),3> _B_bg;
      real_t _rate {};
      real_t(* _profile)(real_t q_normal) = nullptr;

    public:
      auto& set_E_background( apt::array<real_t(*)(real_t,real_t,real_t,real_t t),3> E ) {_E_bg = E; return *this;}
      auto& set_B_background( apt::array<real_t(*)(real_t,real_t,real_t,real_t t),3> B ) {_B_bg = B; return *this;}
      auto& set_damping_rate( real_t rate ) {_rate = rate; return *this;}
      auto& set_damping_profile( real_t(*p)(real_t) ) {_profile = p; return *this;}

      virtual DampingLayer* Clone() const { return new DampingLayer(*this); }

      virtual void operator() ( Field<3>& E, Field<3>& B, JField&,
                                const Grid& grid, const mpi::CartComm&,
                                int timestep, real_t dt) const override {
        auto impose =
          [this, &grid, dt]( auto& F, const auto& F_bg, int comp ) {
            const auto& ofs = F[comp].offset();
            const auto& m = F.mesh();

            for ( int j = range::begin(*this,1); j < range::end(*this,1); ++j ) {
              int li_j = m.linear_index(1, j);
              real_t th = grid[1].absc(j, ofs[1]*0.5);
              for ( int i = range::begin(*this,0); i < range::end(*this,0); ++i ) {
                int li = li_j + i * m.stride(0);
                real_t lnr = grid[0].absc(i, ofs[0] * 0.5);
                real_t lambda = 1.0 - _rate * dt * _profile(lnr);
                real_t f_bg = F_bg[comp](lnr,th,0.0,0.0); // time = 0.0, so damp to the initial condition
                F[comp][li] = ( F[comp][li] - f_bg ) * lambda + f_bg;
              }
            }
          };
        for ( int C = 0; C < 3; ++C )
          impose(E,_E_bg,C);
        for ( int C = 0; C < 3; ++C )
          impose(B,_B_bg,C);
      }
    } fu_damp;
    {
      auto& fu = fu_damp;

      auto profile =
        [](real_t lnr) -> real_t {
          static real_t r_damp_b = std::exp(pic::supergrid[0].upper()) - damping_thickness;
          lnr = ( std::exp(lnr) - r_damp_b ) / damping_thickness;
          return 0.5 * lnr * lnr;
        };

      const int damping_begin_index = supergrid[0].csba(std::log(std::exp(supergrid[0].upper()) - damping_thickness));
      fu.setName("DampingLayer");
      fu[0] = { damping_begin_index, pic::supergrid[0].dim()+myguard };
      fu[1] = { 0, pic::supergrid[1].dim() + 1 };

      fu.set_damping_rate(damping_rate);
      fu.set_damping_profile(profile);
      fu.set_E_background({ E_r_star, E_theta_star, E_phi_star });
      fu.set_B_background({ B_r_star, B_theta_star, B_phi_star });
    }

    struct Axissymmetric : public FieldAction {
    private:
      bool _is_upper_axis = false;

    public:
      auto& is_upper_axis( bool x ) { _is_upper_axis = x; return *this; }
      virtual Axissymmetric* Clone() const { return new Axissymmetric(*this); }
      virtual void operator() ( Field<3>& E, Field<3>& B,
                                JField& , const Grid& grid,
                                const mpi::CartComm&, int , real_t) const override {
        // NOTE Guard cells values are needed when interpolating E and B
        // E_theta, B_r, B_phi are on the axis. All but B_r should be set to zero
        auto assign = []( real_t& v_g, real_t& v_b ) noexcept { v_g = v_b; };
        auto neg_assign = []( real_t& v_g, real_t& v_b ) noexcept {
                            v_g = ( &v_g == &v_b ) ? 0.0_r : - v_b;
                          };
        // MIDWAY in AxisDir
        axissymmetrize(E[0], assign, range::begin(*this), range::end(*this), _is_upper_axis );
        axissymmetrize(E[2], neg_assign, range::begin(*this), range::end(*this), _is_upper_axis );
        axissymmetrize(B[1], neg_assign, range::begin(*this), range::end(*this), _is_upper_axis );

        // INSITU in AxisDir
        axissymmetrize(E[1], neg_assign, range::begin(*this), range::end(*this), _is_upper_axis );
        axissymmetrize(B[0], assign, range::begin(*this), range::end(*this), _is_upper_axis );
        axissymmetrize(B[2], neg_assign, range::begin(*this), range::end(*this), _is_upper_axis );
      }
    } fu_asym_lo, fu_asym_hi;
    {
      fu_asym_lo.setName("AxissymmetrizeEBLower");
      fu_asym_lo[0] = { 0, pic::supergrid[0].dim() };
      fu_asym_lo[1] = { -myguard, 1 }; // NOTE 1 so as to set values right on axis
      fu_asym_lo.is_upper_axis(false);
      fu_asym_lo.require_original_EB(false);

      fu_asym_hi.setName("AxissymmetrizeEBHigher");
      fu_asym_hi[0] = { 0, pic::supergrid[0].dim() };
      fu_asym_hi[1] = { pic::supergrid[1].dim(), pic::supergrid[1].dim() + myguard };
      fu_asym_hi.is_upper_axis(true);
      fu_asym_hi.require_original_EB(false);
    }

    fus.emplace_back(fu.Clone());
    fus.emplace_back(fu_cond.Clone());
    fus.emplace_back(fu_damp.Clone());
    fus.emplace_back(fu_asym_lo.Clone());
    fus.emplace_back(fu_asym_hi.Clone());

    return fus;
  }
}

namespace pic {
  auto set_up_particle_properties() {
    map<Properties> properties;
    {
      properties.insert(species::electron, {1.0,-1.0,"electron","el"});
      properties.insert(species::ion, { 1.0, 1.0, "ion","io"});
    }

    {
      constexpr auto* lorentz = ::particle::force::template lorentz<real_t,Specs,::particle::vParticle>;
      constexpr auto* gravity = ::particle::force::gravity<real_t,Specs,::particle::vParticle>;

      if ( properties.has(species::electron) ) {
        auto sp = species::electron;
        Force force;
        const auto& prop = properties[sp];

        force.add( lorentz, prop.charge_x / prop.mass_x );
        force.add( gravity, gravity_strength );

        force.Register(sp);
      }
      if ( properties.has(species::ion) ) {
        auto sp = species::ion;
        Force force;
        const auto& prop = properties[sp];

        force.add( lorentz, prop.charge_x / prop.mass_x );
        force.add( gravity, gravity_strength );

        force.Register(sp);
      }
    }

    return properties;
  }

}

namespace pic {
  real_t N_atm_floor() {
    return std::exp(1.0) * 2.0 * Omega * mu / wpic2;
  }

  auto set_up_particle_actions() {
    namespace range = apt::range;
    std::vector<std::unique_ptr<PtcAction>> pus;

    struct Updater_with_J_species: PtcAction {
    private:
      PtcUpdater _pu;

    public:
      void set_updater( PtcUpdater&& pu) noexcept { _pu = std::move(pu); }

      Updater_with_J_species* Clone() const override { return new Updater_with_J_species(*this); }

      void operator() ( map<PtcArray>& particles,
                        JField& J,
                        std::vector<Particle>* new_ptc_buf,
                        const map<Properties>& properties,
                        const Field<3>& E,
                        const Field<3>& B,
                        const Grid& grid,
                        const Ensemble* ens,
                        real_t dt, int timestep, util::Rng<real_t>& rng
                        ) override {
        if ( !RTD::data().Jsp || !export_plan.is_do(timestep) ) {
          _pu( particles,J, new_ptc_buf, properties, E, B, grid, ens, dt, timestep, rng );
        } else {
          auto &Jsp = *(RTD::data().Jsp);
          // store J by species separately for data export
          for ( auto sp : particles ) {
            map<PtcArray> ptcs_sp;
            ptcs_sp.insert(sp);
            std::swap( ptcs_sp[sp], particles[sp] );
            Jsp[sp] = JField(J.mesh());
            _pu( ptcs_sp, Jsp[sp], new_ptc_buf, properties, E,B,grid,ens,dt,timestep,rng );
            std::swap( ptcs_sp[sp], particles[sp] );

            for ( int C = 0; C < 3; ++C ) {
              for ( int i = 0; i < J.mesh().linear_size(); ++i )
                J[C][i] += Jsp[sp][C][i];
            }
          }
        }
      }
    } pu;
    {
      PtcUpdater pu0;
      pu0.set_update_q(Metric::geodesic_move<apt::vVec<real_t,3>, apt::vVec<real_t,3>>);

      pu.setName("MainUpdate");
      pu.set_updater(std::move(pu0));
    }

    ::particle::Migrator<DGrid,real_t,Specs,ShapeF,real_j_t> migrate;
    {
      migrate.setName("MigrateParticles");
      migrate.set_supergrid(pic::supergrid);
    }

    struct Atmosphere: PtcAction {
    private:
      Field<1> _count_n;
      Field<1> _count_p;
      int _n = 0; // normal direction
      species _posion = species::ion;
      species _negaon = species::electron;
      real_t _v_th = 0.0;
      real_t _N_atm = 0.0;
      real_t _min_frac = 1e-6; // over fracs larger than this will be injected
      real_t (*_omega_t) ( real_t time ) = nullptr;

    public:
      auto& set_thermal_velocity(real_t v) { _v_th = v; return *this; }
      auto& set_number_in_atmosphere(real_t N) { _N_atm = N; return *this; }
      auto& set_minimal_fraction( real_t x ) { _min_frac = x; return *this; }
      auto& set_normal_direction( int n ) { _n = n; return *this; }
      auto& set_omega_t(real_t (*omega_t) ( real_t )) { _omega_t = omega_t; return *this; }
      auto& set_positive_charge(species sp) { _posion = sp; return *this; }
      auto& set_negative_charge(species sp) { _negaon = sp; return *this; }

      Atmosphere* Clone() const override { return new Atmosphere(*this); }

      void operator() ( map<PtcArray>& particles,
                        JField& J,
                        std::vector<Particle>*,
                        const map<Properties>& properties,
                        const Field<3>& E,
                        const Field<3>& B,
                        const Grid& grid,
                        const Ensemble* ens,
                        real_t dt, int timestep, util::Rng<real_t>& rng
                        ) override {
        if( range::end(*this,_n) <= range::begin(*this,_n) ) return;

        _count_n.resize( {apt::make_range(range::begin(*this),range::end(*this),0)} );
        _count_p.resize( {apt::make_range(range::begin(*this),range::end(*this),0)} );

        apt::array<real_t,DGrid> lb;
        apt::array<real_t,DGrid> ub;
        for ( int i = 0; i < DGrid; ++i ) {
          lb[i] = grid[i].absc( range::begin(*this,i), 0.0 );
          ub[i] = grid[i].absc( range::end(*this,i), 0.0 );
        }
        auto is_in = [&lb,&ub]( const auto& q ) {
                       for ( int i = 0; i < DGrid; ++i ) {
                         if ( q[i] < lb[i] || q[i] >= ub[i] ) return false;
                       }
                       return true;
                     };

        auto f_count
          = [&lb,&ub,&grid,is_in]( auto& count, const auto& ptcs) {
              count.reset();
              for ( const auto& x : ptcs ) {
                if ( !x.is(flag::exist) || !is_in(x.q()) ) continue;
                Index idx;
                for ( int i = 0; i < DGrid; ++i )
                  idx[i] = grid[i].csba(x.q(i)); // NOTE used grid.lower instead of lb, important
                count[0](idx) += x.frac(); // add by fraction
              }
            };

        f_count( _count_n, particles[_negaon] );
        f_count( _count_p, particles[_posion] );

        { // parallelize TODO optimize
          int rank_inj = timestep % ens->size();
          ens->intra.template reduce<true>(mpi::by::SUM, rank_inj, _count_n[0].data().data(), _count_n[0].data().size() );
          ens->intra.template reduce<true>(mpi::by::SUM, rank_inj, _count_p[0].data().data(), _count_p[0].data().size() );
          if ( ens->intra.rank() != rank_inj ) return;
        }

        auto itr_po = std::back_inserter(particles[_posion]);
        auto itr_ne = std::back_inserter(particles[_negaon]);

        for ( const auto& I : apt::Block(range::begin(*this),range::end(*this)) ) {
          auto N_pairs = std::min( _count_n[0](I), _count_p[0](I) );
          Vec3 q{};
          for ( int i = 0; i < DGrid; ++i )
            q[i] = grid[i].absc(I[i], 0.5);

          Vec3 nB {1.0_r, 0.0_r, 0.0_r};

          Vec3 p{};
          p[2] = _omega_t( timestep * dt ) * std::exp(q[0]) * std::sin(q[1]); // corotating

          // replenish
          real_t quota = _N_atm * std::sin(q[1]) - N_pairs;
          while ( quota > _min_frac ) {
            auto q_ptc = q;
            real_t frac = std::min( 1.0_r, quota );
            quota -= 1.0_r;

            for ( int i = 0; i < DGrid; ++i ) {
              if ( _n == i )
                q_ptc[i] += grid[i].delta() * rng.uniform(-0.5, 0.0);
              else
                q_ptc[i] += grid[i].delta() * rng.uniform(-0.5, 0.5);
            }
            auto p_ptc = p;
            p_ptc += nB * rng.gaussian( 0.0, _v_th );
            *(itr_ne++) = Particle( q_ptc, p_ptc, frac, _negaon );
            *(itr_po++) = Particle( std::move(q_ptc), std::move(p_ptc), frac, _posion );
          }
        }
      }
    } atm;
    {
      atm.setName("Atmosphere");
      atm[0] = { star_interior - 1, star_interior };
      atm[1] = { 0, supergrid[1].dim() };

      atm.set_thermal_velocity(v_th).set_number_in_atmosphere(atm_x * N_atm_floor());
      atm.set_positive_charge(species::ion).set_negative_charge(species::electron);
      atm.set_omega_t(omega_spinup).set_normal_direction(0);
    }

    struct Axissymmetric : public PtcAction {
    private:
      bool _is_upper_axis = false;
    public:
      auto& is_upper_axis( bool x ) { _is_upper_axis = x; return *this; }
      virtual Axissymmetric* Clone() const { return new Axissymmetric(*this); }

      virtual void operator() ( map<PtcArray>&, JField& J, std::vector<Particle>*, const map<Properties>&,
                                const Field<3>&, const Field<3>&, const Grid& grid, const Ensemble* ,
                                real_t, int, util::Rng<real_t>&) override {
        auto add_assign =
          []( real_j_t& a, real_j_t& b ) noexcept {
            a += b;
            b = a;
          };

        auto sub_assign =
          []( real_j_t& a, real_j_t& b ) noexcept {
            a -= b;
            b = -a;
          };
        // MIDWAY in AxisDir
        axissymmetrize(J[0], add_assign, range::begin(*this),range::end(*this),_is_upper_axis );
        axissymmetrize(J[2], sub_assign, range::begin(*this),range::end(*this),_is_upper_axis );
        // INSITU in AxisDir
        axissymmetrize(J[1], sub_assign, range::begin(*this),range::end(*this),_is_upper_axis );
      }
    } asym_lo, asym_hi;
    {
      asym_lo.setName("AxissymmetrizeJLower");
      asym_lo[0] = { 0, supergrid[0].dim() };
      asym_lo[1] = { -myguard, 0 }; // NOTE 0, cannot be 1 because it's add/sub_assign
      asym_lo.is_upper_axis(false);

      asym_hi.setName("AxissymmetrizeJHigher");
      asym_hi[0] = { 0, supergrid[0].dim() };
      asym_hi[1] = { supergrid[1].dim(), supergrid[1].dim() + myguard };
      asym_hi.is_upper_axis(true);
    }

    struct Escaping : public PtcAction {
    private:
      int _n = 0; // normal direction

    public:
      virtual Escaping* Clone() const { return new Escaping(*this); }

      virtual void operator() ( map<PtcArray>& particles, JField& J, std::vector<Particle>* new_ptc_buf, const map<Properties>&,
                                const Field<3>&, const Field<3>&, const Grid& grid, const Ensemble* ,
                                real_t dt, int timestep, util::Rng<real_t>& rng) override {
        if ( apt::range::is_empty(*this) ) return;

        for ( auto sp : particles ) {
          for ( auto ptc : particles[sp] ) { // TODOL semantics
            // NOTE: this implementation assumes that there is never particles traveling backwards. If not, the flag::ignore_force must be unset when such particles leave the layer, and one must be careful that the ignore_force delimiter must not coincide with the lower bound of a patch.
            if ( ptc.q(_n) >= grid[_n].absc(range::begin(*this,_n)) and ptc.q(_n) < grid[_n].absc(range::end(*this,_n)) ) {
              if ( not ptc.is(flag::ignore_force) ) {
                // turn momentum to radial
                auto p = apt::abs(ptc.p());
                ptc.p(1) = 0;
                ptc.p(2) = 0;
                ptc.p(_n) = p;
                ptc.set(flag::ignore_force);
              }
            }
          }
        }
      }
    } escape;
    {
      escape.setName("Escaping");
      const int escape_begin_index = supergrid[0].csba(std::log(std::exp(supergrid[0].upper()) - 0.5*damping_thickness)); // at halfway so as to reduce reflection of artifacts
      escape[0] = { escape_begin_index, supergrid[0].dim() + myguard };
      escape[1] = { 0, supergrid[1].dim() };
    }

    pus.emplace_back(escape.Clone());
    pus.emplace_back(pu.Clone());
    pus.emplace_back(atm.Clone());
    pus.emplace_back(asym_lo.Clone());
    pus.emplace_back(asym_hi.Clone());
    // FIXME migrate need more memory check
    pus.emplace_back(migrate.Clone()); // After this line, particles are all within borders.

    return pus;
  }
}

namespace pic {
  auto set_up_initial_conditions() {
    // local class in a function
    struct InitialCondition : public apt::ActionBase<DGrid> {
      InitialCondition* Clone() const override { return new InitialCondition(*this); }

      void operator() ( const Grid& grid,
                        Field<3>& ,
                        Field<3>& B,
                        JField& ,
                        map<PtcArray>&
                        ) const {
        for ( const auto& I : apt::Block(apt::range::begin(*this),apt::range::end(*this)) ) {
          B[0](I) = B_r_star( grid[0].absc(I[0], 0.5 * B[0].offset()[0]), grid[1].absc(I[1], 0.5 * B[0].offset()[1]), 0, 0 );
        }
      }

    } ic;
    ic[0] = { 0, supergrid[0].dim() + 1 };
    ic[1] = { 0, supergrid[1].dim() + 1 }; // NOTE +1 to include upper boundary

    return ic;
  }

  auto set_up_post_resume_actions() {
    struct PostResume : public apt::ActionBase<DGrid> {
      PostResume *Clone() const override {return new PostResume(*this);}

      void
      operator()(const Grid &grid, Field<3> &E, Field<3> &B, JField &J,
                 map<PtcArray> &particles,
                 const std::optional<Ensemble> & ens_opt,
                 int resumed_timestep,
                 std::string this_run_dir) const {}
    } pr;
    pr[0] = {0, supergrid[0].dim() + myguard};
    pr[1] = {0, supergrid[1].dim() + 1}; // NOTE +1 to include upper boundary

    return pr;
  }
}

namespace pic {
  constexpr bool is_collinear_mesh = false; // FIXME this is an ad hoc fix

  void export_prior_hook( const map<PtcArray>& particles, const map<Properties>& properties,
                          const Field<3>& E, const Field<3>& B, const JField& J,  const Grid& grid, const std::optional<mpi::CartComm>& cart_opt, const Ensemble& ens,
                          real_t dt, int timestep ) {
    if (RTD::data().Jsp) {
      auto& Jsp = *RTD::data().Jsp;
      for ( auto s : Jsp ) {
        for ( int i = 0; i < 3; ++i )
          ens.reduce_to_chief( mpi::by::SUM, Jsp[s][i].data().data(), Jsp[s][i].data().size() );
      }
      if ( cart_opt ) {
        for ( auto s : Jsp ) {
          field::merge_sync_guard_cells( Jsp[s], *cart_opt );
        }
      }
    }
  }
}

namespace pic {
  using RDS = real_export_t;
  using IOField = ::field::Field<RDS,3,DGrid>;
  using IOGrid = ::apt::Grid<RDS,DGrid>;

  constexpr auto I2std ( const Index& I ) {
    apt::array<real_t, DGrid> res;
    for ( int i = 0; i < DGrid; ++i )
      res[i] = I[i] + 0.5; // interpolate to MIDWAY
    return res;
  }

  template <int F>
  apt::array<real_t,3> field_self ( Index I, const Grid& grid, const Field<3>& E,
                                    const Field<3>& B, const JField& J ) {
    if constexpr ( F == 0 ) {
        return msh::interpolate( E, I2std(I), ShapeF() );
      } else if ( F == 1 ) {
      return msh::interpolate( B, I2std(I), ShapeF() );
    } else if ( F == 2 ) {
      return msh::interpolate( J, I2std(I), ShapeF() );
    } else {
      static_assert(F < 3);
    }
  }

  constexpr int POW( int B, int E ) {
    if ( E == 0 ) return 1;
    else return B * POW(B,E-1);
  }

  void average_when_downsampled ( IOField& fds, const IOGrid& grid, int num_comps, const mpi::CartComm& ) {
    int downsample_ratio =
        static_cast<int>(grid[0].delta() / pic::supergrid[0].delta() + 0.5);
    int factor = POW(downsample_ratio, pic::DGrid);
    for ( int i = 0; i < num_comps; ++i ) {
      for ( auto& x : fds[i].data() ) x /= factor;
    }
  }

  void average_and_divide_flux_by_area ( IOField& fds, const IOGrid& grid, int num_comps, const mpi::CartComm& cart ) {
    average_when_downsampled(fds, grid, num_comps, cart);

    using Metric = metric::LogSpherical<RDS>;
    // define a function pointer.
    RDS (*hh_func)(RDS,RDS,RDS) = nullptr;
    apt::array<RDS,3> q {};

    for ( int comp = 0; comp < num_comps; ++comp ) {
      const auto& ofs = fds[comp].offset();
      switch(comp) {
      case 0: hh_func = Metric::template hh<0>; break;
      case 1: hh_func = Metric::template hh<1>; break;
      case 2: hh_func = Metric::template hh<2>; break;
      }

      for ( const auto& I : apt::Block(apt::range::begin(fds.mesh().range()), apt::range::end(fds.mesh().range())) ) {
        for ( int i = 0; i < DGrid; ++i ) q[i] = grid[i].absc(I[i], 0.5 * ofs[i]);
        auto hh = hh_func(q[0], q[1], q[2]);
        if ( std::abs(hh) > 1e-12 )
          fds[comp](I) /= hh;
        else
          fds[comp](I) = 0.0;
      }
    }
  }

  template < particle::species SP >
  apt::array<real_t,3> J_by_species ( Index I, const Grid& grid, const Field<3>& ,
                                      const Field<3>& , const JField& ) {
    // NOTE due to interpolation, the exported Jsp's don't sum up to J.
    const auto& Jsp = (*RTD::data().Jsp)[SP];
    return msh::interpolate(Jsp, I2std(I), ShapeF());
  }
}

namespace pic {
  apt::array<real_t,3> ptc_num ( const Properties& prop, const typename PtcArray::const_particle_type& ptc ) {
    return { 1.0, 0.0, 0.0 };
  }

  // NOTE energy should factor in mass
  apt::array<real_t,3> ptc_energy ( const Properties& prop, const typename PtcArray::const_particle_type& ptc ) {
    if (prop.mass_x > 0.01_r)
      return { prop.mass_x * std::sqrt( 1.0_r + apt::sqabs(ptc.p()) ), 0.0, 0.0 };
    else
      return { apt::abs(ptc.p()), 0.0, 0.0};
  }

  apt::array<real_t,3> ptc_momentum ( const Properties& prop, const typename PtcArray::const_particle_type& ptc ) {
    if (prop.mass_x > 0.01_r)
      return { prop.mass_x*ptc.p(0), prop.mass_x*ptc.p(1), prop.mass_x*ptc.p(2) };
    else
      return { ptc.p(0), ptc.p(1), ptc.p(2) };
  }
}

namespace pic {

  using DataExporter_t = io::DataExporter<real_export_t, DGrid, real_t, particle::Specs, real_j_t>;

  std::vector<DataExporter_t> set_up_data_exporters() {
    auto add_common_exportees = [] ( auto& exporter ) {
      using FA = ::io::FexpTbyFunction<real_export_t, DGrid, real_t, real_j_t>;

      exporter
        .add_exportee( new FA ( "E", 3, field_self<0>, average_when_downsampled) )
        .add_exportee( new FA ( "B", 3, field_self<1>, average_when_downsampled) )
        .add_exportee( new FA ( "J", 3, field_self<2>, average_and_divide_flux_by_area) );

      if ( RTD::data().Jsp ) {
        using namespace particle;
        for (auto sp : *(RTD::data().Jsp)) {
          switch(sp) {
          case species::electron :
            exporter.add_exportee(new FA(
                                              "Je", 3, J_by_species<species::electron>,
                                              average_and_divide_flux_by_area));
            break;
          case species::positron :
            exporter.add_exportee(new FA(
                                              "Jp", 3, J_by_species<species::positron>,
                                              average_and_divide_flux_by_area));
            break;
          case species::ion :
            exporter.add_exportee(new FA(
                                              "Ji", 3, J_by_species<species::ion>,
                                              average_and_divide_flux_by_area));
            break;
          default: ;
          }
        }
      }
      using PA = ::io::PexpTbyFunction<real_export_t, DGrid, real_t, Specs>;
      exporter
        .add_exportee(new PA("Num", 1, ptc_num, nullptr))
        .add_exportee(new PA("E", 1, ptc_energy, nullptr))
        .add_exportee(new PA("P", 3, ptc_momentum, nullptr));
    };


    // Note that DataExporter_t doesn't have copy constructor
    std::vector<DataExporter_t> res;
    {
      auto& main_exporter = res.emplace_back();
      main_exporter.set_is_collinear_mesh(false)
        .set_downsample_ratio(export_plan.downsample_ratio)
        .set_data_dir("")
        .set_range(
            {{{0, supergrid[0].dim() + myguard}, {0, supergrid[1].dim() + 1}}});
      add_common_exportees(main_exporter);
    }

    return res;
  }
}

namespace pic {
  void export_post_hook() {
    if ( RTD::data().Jsp ) {
      // clear Jsp to save some space
      for (auto sp : *(RTD::data().Jsp))
        (*RTD::data().Jsp)[sp] = {};
    }
  }
}

#include <sstream>
#include "apt/print.hpp"

namespace pic {
  std::string proofread(std::string indent) {
    std::ostringstream o;
    real_t gamma_0 = Omega * Omega * mu;
    o << indent << "Np=" << apt::fmt("%.1f", 2 * Omega * mu / wpic2 ) << std::endl;
    o << indent << "(w_pic dt)^2 = " << apt::fmt("%.4f", wpic2 * dt * dt ) << std::endl;
    o << indent << "re=" << apt::fmt("%.4f", r_e() ) << std::endl;

    {
      using namespace particle;
      o << indent << "ATM: N_atm_floor=" << apt::fmt("%.1f", N_atm_floor());
      o << ", atm_x=" << apt::fmt("%.1f", atm_x);
      o << ", v_th=" << apt::fmt("%.2f", v_th)
        << ", g=" << apt::fmt("%.2f", gravity_strength) << std::endl;
    }

    return o.str();
  }
}

#endif
