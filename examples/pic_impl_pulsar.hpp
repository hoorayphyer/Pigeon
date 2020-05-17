#ifndef _PIC_IMPL_HPP_
#define _PIC_IMPL_HPP_

#include "../examples/pic_prior_impl.hpp"

#include "apt/numeric.hpp"

#include "particle/annihilation.hpp"

#include "pic/plans.hpp"
#include "pic/forces/gravity.hpp"
#include "pic/forces/landau0.hpp"

#include "metric/log_spherical.hpp"

#include "field/old_field_solver/updater.hpp" // FIXME

#include "io/exportee_by_function.hpp"

namespace pic {
  using Metric = metric::LogSpherical<real_t>;

  inline constexpr apt::array<int,DGrid> dims = { 1, 1 };
  inline constexpr apt::array<bool,DGrid> periodic = {false,false};

  constexpr Grid supergrid
  = {{ { 0.0, std::log(30.0), 64 }, { 0.0, 180.0_deg, 64 } }};

  real_t Omega = 1.0 / 6.0;
  real_t dt;
  real_t mu;
  real_t wpic2;

  real_t gamma_fd;
  real_t E_ph;
  real_t gravity_strength;
  real_t landau0_B_thr;
  real_t magnetic_convert_B_thr;
  std::array<real_t,2> mfp;
  real_t atm_x;
  real_t v_th;

  int damping_layer;
  real_t damping_rate;
  real_t spinup_time;

  std::string project_name;
  std::string datadir_prefix;
  int total_timesteps;

  Plan annih_plan {};

  void load_configuration(const ConfFile_t& conf) {
    safe_set(dt, conf["dt"]);
    safe_set(mu, conf["gamma0"]); mu /= std::pow(Omega,2.0);
    safe_set(wpic2, conf["Np"]); wpic2 = 2 * Omega * mu / wpic2;

    safe_set(gamma_fd,conf["pairs"]["gamma_fd"]);
    safe_set(E_ph, conf["pairs"]["E_ph"]);
    safe_set(gravity_strength, conf["forces"]["gravity"]);
    safe_set(landau0_B_thr, conf["forces"]["landau0_ratio"]); landau0_B_thr *= pic::mu;

    safe_set(magnetic_convert_B_thr, conf["pairs"]["photon"]["magnetic_convert_ratio"]); magnetic_convert_B_thr *= pic::mu;
    safe_set(mfp[0], conf["pairs"]["photon"]["mfp"][0]);
    safe_set(mfp[1], conf["pairs"]["photon"]["mfp"][1]);

    safe_set(v_th, conf["atmosphere"]["v_th"]);
    safe_set(atm_x, conf["atmosphere"]["multiplicity"]);

    safe_set(damping_layer, conf["damping"]["layer"]);
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

    save_tracing_plan.on = conf["plans"]["tracing"]["on"].value_or(false);
    save_tracing_plan.start = conf["plans"]["tracing"]["start"].value_or(0);
    safe_set(save_tracing_plan.interval, conf["plans"]["tracing"]["interval"]);
    save_tracing_plan.num_files = conf["plans"]["tracing"]["num_files"].value_or(1);

    annih_plan.on = conf["plans"]["annihilation"]["on"].value_or(false);
    annih_plan.start = conf["plans"]["annihilation"]["start"].value_or(0);
    safe_set(annih_plan.interval, conf["plans"]["annihilation"]["interval"]);
  }

  real_t r_e() {
    real_t res = wpic2 / ( 4.0 * std::acos(-1.0l));
    res *= apt::dV(supergrid);
    return res;
  }
}

namespace pic {
  constexpr int star_interior = 5;

  constexpr int myguard = std::max(1, ( pic::ShapeF::support() + 3 ) / 2 ); // NOTE minimum number of guards of J on one side is ( supp + 3 ) / 2

  real_t omega_spinup ( real_t time ) noexcept {
    return std::min<real_t>( time / spinup_time, 1.0_r ) * Omega;
  }

  real_t B_r_star( real_t lnr, real_t theta, real_t , real_t time ) noexcept {
    return pic::mu * 2.0_r * std::cos(theta) * std::exp(-3.0_r * lnr);
  }
  real_t B_theta_star( real_t lnr, real_t theta, real_t , real_t time ) noexcept {
    return pic::mu * std::sin(theta) * std::exp(-3.0_r * lnr);
  }
  constexpr real_t B_phi_star( real_t lnr, real_t theta, real_t , real_t time ) noexcept {
    return 0;
  }

  real_t E_r_star( real_t lnr, real_t theta, real_t , real_t time ) noexcept {
    return pic::mu * omega_spinup(time) * std::exp(- 2.0_r * lnr) * std::sin( theta ) * std::sin(theta);
  }
  real_t E_theta_star( real_t lnr, real_t theta, real_t , real_t time ) noexcept {
    return - pic::mu * omega_spinup(time) * std::exp(- 2.0_r * lnr ) * std::sin( 2.0_r*theta );
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
  // TODOL annihilation will affect deposition // NOTE one can deposit in the end
}

namespace pic {
  struct RTD { // runtime data
  public:
    static RTD& data() {
      static RTD r;
      return r;
    }

    map<real_t> N_scat {};
    Field<1> pc_counter {};
    real_t pc_cumulative_time {};

    map<JField> Jsp {}; // current by species
    bool is_export_Jsp = false;

    Field<1> skin_depth {};

    void init( const map<Properties>& properties, const Grid& localgrid ) {
      is_export_Jsp = false;
      for ( auto sp : properties )
        N_scat.insert( sp, 0 );
      if ( is_export_Jsp ) {
        for ( auto sp : properties )
          Jsp.insert( sp, {} );
      }

      Index bulk_dims;
      for ( int i = 0; i < DGrid; ++i ) bulk_dims[i] = localgrid[i].dim();
      auto range = apt::make_range({}, bulk_dims, myguard); // FIXME range with no guard gives memory error. Interpolation in export needs them.
      pc_counter = {range};
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

    ::field::OldSolve<real_t,DGrid,real_j_t> fu;
    {
      const int guard = fu.guard();
      fu.setName("OldSolve");
      fu[0] = { 0, pic::supergrid[0].dim(), guard };
      fu[1] = { 0, pic::supergrid[1].dim(), guard };

      fu.set_fourpi( 4.0 * std::acos(-1.0l) * pic::r_e() );
      fu.set_mu(pic::mu);
      fu.set_magnetic_pole(2);
      fu.set_damping_rate(damping_rate);
      fu.set_surface_indent(5);
      fu.set_damp_indent(damping_layer);
      fu.set_omega_t(omega_spinup);
    }

    struct Axissymmetric : public FieldAction {
    private:
      bool _is_upper_axis = false;

    public:
      auto& is_upper_axis( bool x ) { _is_upper_axis = x; return *this; }
      virtual Axissymmetric* Clone() const { return new Axissymmetric(*this); }
      virtual void operator() ( Field<3>& E, Field<3>& B,
                                const JField& , const Grid& grid,
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
    fus.emplace_back(fu_asym_lo.Clone());
    fus.emplace_back(fu_asym_hi.Clone());

    return fus;
  }
}

namespace pic {
  constexpr real_t gamma_off = 15.0;
  constexpr real_t Ndot_fd = 0.25;

  auto set_up_particle_properties() {
    map<Properties> properties;
    {
      properties.insert(species::electron, {1.0,-1.0,"electron","el"});
      properties.insert(species::positron, {1.0,1.0,"positron","po"});
      properties.insert(species::ion, { 5.0, 1.0, "ion","io"});
      properties.insert(species::photon, { 0, 0, "photon","ph" });
    }

    {
      constexpr auto* lorentz = ::particle::force::template lorentz<real_t,Specs,::particle::vParticle>;
      constexpr auto* landau0 = ::particle::force::landau0<real_t,Specs,::particle::vParticle>;
      constexpr auto* gravity = ::particle::force::gravity<real_t,Specs,::particle::vParticle>;

      if ( properties.has(species::electron) ) {
        auto sp = species::electron;
        Force force;
        const auto& prop = properties[sp];

        force.add( lorentz, prop.charge_x / prop.mass_x );
        force.add( gravity, gravity_strength );
        force.add( landau0, landau0_B_thr );

        force.Register(sp);
      }
      if ( properties.has(species::positron) ) {
        auto sp = species::positron;
        Force force;
        const auto& prop = properties[sp];

        force.add( lorentz, prop.charge_x / prop.mass_x );
        force.add( gravity, gravity_strength );
        force.add( landau0, landau0_B_thr );

        force.Register(sp);
      }
      if ( properties.has(species::ion) ) {
        auto sp = species::ion;
        Force force;
        const auto& prop = properties[sp];

        force.add( lorentz, prop.charge_x / prop.mass_x );
        force.add( gravity, gravity_strength );
        // force.add( landau0, landau0_B_thr );

        force.Register(sp);
      }
    }

    using Ptc_t = typename PtcArray::particle_type;
    namespace scat = ::particle::scat;
    {
      ::particle::Scat<real_t,Specs> ep_scat;

      ep_scat.eligs.push_back([](const Ptc_t& ptc){ return ptc.q(0) < std::log(9.0_r); });

      scat::CurvatureRadiation<real_t,Specs>::gamma_fd = gamma_fd;
      scat::CurvatureRadiation<real_t,Specs>::gamma_off = gamma_off;
      scat::CurvatureRadiation<real_t,Specs>::Ndot_fd = Ndot_fd;
      scat::CurvatureRadiation<real_t,Specs>::E_ph = E_ph;
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
      ::particle::Scat<real_t,Specs> photon_scat;
      // Photons are free to roam across all domain. They may produce pairs outside light cylinder
      photon_scat.eligs.push_back([](const Ptc_t& ptc) { return true; });
      scat::MagneticConvert<real_t,Specs>::B_thr = magnetic_convert_B_thr;
      scat::MagneticConvert<real_t,Specs>::mfp = mfp[0];
      photon_scat.channels.push_back( scat::MagneticConvert<real_t,Specs>::test );

      scat::TwoPhotonCollide<real_t,Specs>::mfp = mfp[1];
      photon_scat.channels.push_back( scat::TwoPhotonCollide<real_t,Specs>::test );

      photon_scat.impl = scat::PhotonPairProduction<real_t,Specs>;

      photon_scat.Register( species::photon );
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
        if ( !RTD::data().is_export_Jsp || !export_plan.is_do(timestep) ) {
          _pu( particles,J, new_ptc_buf, properties, E, B, grid, ens, dt, timestep, rng );
        } else {
          auto& Jsp = RTD::data().Jsp;
          // store J by species separately for data export
          for ( auto sp : particles ) {
            map<PtcArray> ptcs_sp;
            ptcs_sp.insert(sp);
            std::swap( ptcs_sp[sp], particles[sp] );
            Jsp[sp] = J;
            Jsp[sp].reset();
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

          Vec3 nB {};
          { // make nB centered in the cell
            const auto& m = B.mesh();
            auto li = m.linear_index(I);
            if constexpr (DGrid == 2) {
                nB[0] = 0.5_r * ( B[0][li] + B[0][li + m.stride(1)] );
                nB[1] = 0.5_r * ( B[1][li] + B[1][li + m.stride(0)] );
                nB[2] = 0.25_r * ( B[2][li] + B[2][li + m.stride(0)] + B[2][li + m.stride(1)] + B[2][li + m.stride(0) + m.stride(1)] );
              } else if (DGrid == 3){
              nB[0] = 0.25_r * ( B[0][li] + B[0][li + m.stride(1)] + B[0][li + m.stride(2)] + B[0][li + m.stride(1) + m.stride(2)] );
              nB[1] = 0.25_r * ( B[1][li] + B[1][li + m.stride(2)] + B[1][li + m.stride(0)] + B[1][li + m.stride(2) + m.stride(0)] );
              nB[2] = 0.25_r * ( B[2][li] + B[2][li + m.stride(0)] + B[2][li + m.stride(1)] + B[2][li + m.stride(0) + m.stride(1)] );
            }
            if ( apt::abs(nB) == 0.0_r ) nB = {1.0_r, 0.0_r, 0.0_r}; // use radial direction as default
            else nB /= apt::abs(nB);
          }

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
                q_ptc[i] += grid[i].delta() * rng.uniform(-0.5, 0.0); // TODO move this out. This only affects when it is lower
              else
                q_ptc[i] += grid[i].delta() * rng.uniform(-0.5, 0.5);
            }
            auto p_ptc = p;
            p_ptc += nB * rng.gaussian( 0.0, _v_th );
            *(itr_ne++) = Particle( q_ptc, p_ptc, frac, _negaon, ::particle::birthplace(ens->label()) );
            *(itr_po++) = Particle( std::move(q_ptc), std::move(p_ptc), frac, _posion, ::particle::birthplace(ens->label()) );
          }
        }
      }
    } atm;
    {
      atm.setName("Atmosphere");
      atm[0] = { star_interior, star_interior + 1 };
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
      asym_lo[1] = { -myguard, 1 }; // NOTE +1 so as to set values right on axis
      asym_lo.is_upper_axis(false);

      asym_hi.setName("AxissymmetrizeJHigher");
      asym_hi[0] = { 0 , supergrid[0].dim() };
      asym_hi[1] = { supergrid[1].dim(), supergrid[1].dim() + myguard };
      asym_hi.is_upper_axis(true);
    }

    struct NewPtcAnalyzer : public PtcAction {
      virtual NewPtcAnalyzer* Clone() const { return new NewPtcAnalyzer(*this); }

      virtual void operator() ( map<PtcArray>& particles, JField& J, std::vector<Particle>* new_ptc_buf, const map<Properties>&,
                                const Field<3>&, const Field<3>&, const Grid& grid, const Ensemble* ,
                                real_t dt, int timestep, util::Rng<real_t>& rng) override {
        // Put particles where they belong after scattering
        assert(new_ptc_buf != nullptr);

        for ( auto& ptc : *new_ptc_buf ) {
          if ( not ptc.is(flag::exist) ) continue;
          const auto this_sp = ptc.template get<species>();

          if ( ptc.is(flag::secondary) ) {
            // log scattering events
            RTD::data().N_scat[this_sp] += ptc.frac();

            // log pair creation events
            if ( species::electron == this_sp ) {
              Index I; // domain index, not the global index
              for ( int j = 0; j < DGrid; ++j )
                I[j] = grid[j].csba( ptc.q(j) );
              RTD::data().pc_counter[0](I) += ptc.frac();
            }
          }

          particles[this_sp].push_back(std::move(ptc));
        }
        new_ptc_buf->resize(0);
        RTD::data().pc_cumulative_time += dt;

      }
    } analyzer;
    {
      analyzer.setName("NewPtcAnalysis");
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
      escape[0] = { supergrid[0].dim() - pic::damping_layer, supergrid[0].dim() + myguard };
      escape[1] = { 0, supergrid[1].dim() };
    }

    class Annihilator : public PtcAction {
    private:
      // policy returns number of pairs to be annihilated
      real_t(*_policy)(real_t num_electron_in_a_cell, real_t num_positron_in_a_cell) = nullptr;

    public:
      Annihilator* Clone() const override { return new Annihilator(*this); }

      Annihilator& set_policy( real_t(*policy)(real_t,real_t) ) noexcept { _policy = policy; return *this; }

      void operator() ( map<PtcArray>& particles, JField& J,
                        std::vector<Particle>* ,
                        const map<Properties>& properties,
                        const Field<3>&, const Field<3>&,
                        const Grid& grid, const Ensemble* ens,
                        real_t dt, int timestep, util::Rng<real_t>& ) override {
        if ( apt::range::is_empty(*this) ) return;

        if ( annih_plan.is_do(timestep) ) {
          assert( ens != nullptr );
          assert( _policy != nullptr );

          apt::array<apt::array<real_t,2>,DGrid> bounds;
          for ( int i = 0; i < DGrid; ++i ) {
            // NOTE guard cells or margins of range not included.
            bounds[i][0] = grid[i].absc( range::begin(*this,i) );
            bounds[i][1] = grid[i].absc( range::end(*this,i) );
          }

          particle::annihilate(particles[species::electron],particles[species::positron],J,
                               properties[species::electron].charge_x,properties[species::positron].charge_x,
                               grid, ens->intra, dt, ShapeF(), _policy, bounds);
        }
      }

    } annih;

    {
      annih.setName("Annihilation");
      annih[0] = { supergrid[0].csba(std::log(9.0)), supergrid[0].dim() + myguard };
      int i_equator = supergrid[1].dim() / 2 + 1; // pick the cell whose lb is at equator.
      int half_width = ( 15.0_deg ) / supergrid[1].delta();
      annih[1] = { i_equator - half_width, i_equator + half_width };

      auto policy =
        []( real_t num_e, real_t num_p ) noexcept {
          return std::min(num_e, num_p) / 2;
        };

      annih.set_policy(policy);
    }

    class NoPhotonZone : public PtcAction {
    private:
      int _interval = 100;

    public:
      NoPhotonZone* Clone() const override { return new NoPhotonZone(*this); }
      NoPhotonZone& set_interval( int interval ) noexcept { _interval = interval; return *this; }

      void operator() ( map<PtcArray>& particles, JField& ,
                        std::vector<Particle>* ,
                        const map<Properties>& ,
                        const Field<3>&, const Field<3>&,
                        const Grid& grid, const Ensemble* ,
                        real_t , int timestep, util::Rng<real_t>& ) override {
        if ( (timestep % _interval != 0) or apt::range::is_empty(*this) ) return;

        apt::array< apt::array<real_t,2>,DGrid> bounds;
        for ( int i = 0; i < DGrid; ++i ) {
          bounds[i][0] = grid[i].absc(range::begin(*this,i));
          bounds[i][1] = grid[i].absc(range::end(*this,i));
        }

        auto in_zone =
          [&bounds] ( const auto& q ) noexcept {
            for ( int i = 0; i < DGrid; ++i ) {
              if ( q[i] < bounds[i][0] or q[i] >= bounds[i][1] )
                return false;
            }
            return true;
          };

        for ( auto ptc : particles[species::photon] ) { // TODOL semantics
          if ( ptc.is(flag::exist) and in_zone(ptc.q()) ) {
            ptc.reset(flag::exist);
          }
        }
      }

    } no_ph;

    {
      no_ph.setName("NoPhotonZone");
      no_ph[0] = { supergrid[0].csba(std::log(9.0)), supergrid[0].dim() + myguard };
      int i_equator = supergrid[1].dim() / 2 + 1; // pick the cell whose lb is at equator.
      int half_width = ( 15.0_deg ) / supergrid[1].delta();
      no_ph[1] = { i_equator - half_width, i_equator + half_width };

      no_ph.set_interval(100);
    }

    struct Tracer : public PtcAction {
    private:
      real_t _prob = 1.0_r;
      std::vector<::particle::species> _sps;

      bool _is_check_within_range = true;

      bool _is_within_range(const PtcArray::particle_type::vec_type& q,
                            const apt::array<apt::array<real_t,2>,DGrid>& bds ) {
        for ( int i = 0; i < DGrid; ++i ) {
          if ( q[i] < bds[i][0] or q[i] >= bds[i][1] ) return false;
        }
        return true;
      }

      using FCond_t = bool(*)(const PtcArray::particle_type& ptc);
      FCond_t _conditional = nullptr;

      Plan _plan{};

      using FMark_t = void(*)(PtcArray::particle_type& ptc);
      FMark_t _marker = nullptr;

    public:
      Tracer* Clone() const override {return new auto(*this);}

      auto& set_probability( real_t prob ) noexcept { _prob = prob; return *this;}
      auto& set_marker( FMark_t f ) noexcept { _marker = f; return *this;}
      auto& set_species(const std::vector<::particle::species>& sps) noexcept {
        _sps = sps; return *this;
      }
      auto& set_is_check_within_range( bool a ) noexcept {_is_check_within_range=a; return *this;}
      auto& set_conditional( FCond_t cond) noexcept {_conditional = cond; return *this;}
      auto& set_plan( const Plan& p ) noexcept { _plan = p; return *this; }

      void operator() ( map<PtcArray>& particles, JField& ,
                        std::vector<Particle>* ,
                        const map<Properties>& ,
                        const Field<3>&, const Field<3>&,
                        const Grid& grid, const Ensemble* ,
                        real_t , int timestep, util::Rng<real_t>& rng) override {
        if ( !_plan.is_do(timestep) or apt::range::is_empty(*this) or !_marker ) return;

        apt::array< apt::array<real_t,2>, DGrid > bds;
        for ( int i = 0; i < DGrid; ++i ) {
          bds[i][0] = grid[i].absc( apt::range::begin(*this,i) );
          bds[i][1] = grid[i].absc( apt::range::end(*this,i) );
        }

        for ( auto sp : _sps ) {
          for ( auto ptc : particles[sp] ) { // TODOL semantics
            if ( !ptc.is(flag::exist)
                 or (_is_check_within_range and !_is_within_range(ptc.q(), bds) )
                 or ( _conditional and !_conditional(ptc) )
                 or ( _prob < 1.0_r and rng.uniform() > _prob )
                 ) continue;
            _marker(ptc);
          }
        }
      }
    };

    Tracer sep_ftp;
    {
      auto& tr = sep_ftp;
      tr.setName("Separatrix Footpoint Tracer");
      tr[0] = { supergrid[0].csba(std::log(1.20)), supergrid[0].csba(std::log(1.40)) };
      tr[1] = { supergrid[1].csba(0.595), supergrid[1].csba(0.638) };
    }

    Tracer ycloud;
    {
      auto& tr = ycloud;
      tr.setName("Y Cloud Tracer");
      tr[0] = { supergrid[0].csba(std::log(4.6)), supergrid[0].csba(std::log(5.2)) };
      tr[1] = { supergrid[1].csba( 90.0_deg - 0.04), supergrid[1].csba(90.0_deg + 0.04) };
    }

    Tracer two_small_plasmoids;
    {
      auto& tr = two_small_plasmoids;
      tr.setName("Two Small Plasmoids Tracer");
      tr[0] = { supergrid[0].csba(std::log(6.6)), supergrid[0].csba(std::log(7.8)) };
      tr[1] = { supergrid[1].csba( 90.0_deg - 0.028), supergrid[1].csba(90.0_deg + 0.028) };
    }

    Tracer one_big_plasmoid;
    {
      auto& tr = one_big_plasmoid;
      tr.setName("One Big Plasmoid Tracer");
      tr[0] = { supergrid[0].csba(std::log(9.6)), supergrid[0].csba(std::log(10.2)) };
      tr[1] = { supergrid[1].csba( 90.0_deg - 0.06), supergrid[1].csba(90.0_deg + 0.06) };
    }

    Tracer black_region;
    {
      auto& tr = black_region;
      tr.setName("Black Region Tracer");
      tr[0] = { supergrid[0].csba(std::log(9.6)), supergrid[0].csba(std::log(10.2)) };
      tr[1] = { supergrid[1].csba(0.582), supergrid[1].csba(0.810) };

      tr.set_is_check_within_range(false);
      tr.set_conditional
        ([](const auto& ptc){
           constexpr apt::array<apt::array<real_t,2>, 4> vs
             = {{ {{0.975,1.478}}, {{1.104,1.280}}, {{1.576,1.500}}, {{1.436, 1.750}} }}; // vertices in x,z cartesian coordinates
           // if a point is inside this polygon, the cross product should all have positive sign
           real_t x = std::exp(ptc.q(0));
           real_t z = x * std::cos(ptc.q(1));
           x = std::sqrt(x*x - z*z);
           for ( int i = 0; i < 4; ++i ) {
             if ( (vs[i][0]-x)*(vs[(i+1)%4][1]-vs[i][1])
                  - (vs[i][1]-z)*(vs[(i+1)%4][0]-vs[i][0]) < 0 )
               return false;
           }
           return true;
         });
    }

    Tracer backflow;
    {
      auto& tr = backflow;
      tr.setName("Backflow Tracer");
      tr[0] = { supergrid[0].csba(std::log(6.0)), supergrid[0].csba(std::log(30)) };
      tr[1] = { supergrid[1].csba(90.0_deg - 0.139), supergrid[1].csba(90.0_deg + 0.139) };
    }
    Tracer grand_tot;
    {
      auto& tr = grand_tot;
      tr.setName("Grand Total");
      tr[0] = { 0, supergrid[0].dim() };
      tr[1] = { 0, supergrid[1].dim() + 1 };
    }


    Plan p_onetime;
    p_onetime.on = true; p_onetime.start = save_tracing_plan.start, p_onetime.interval = 100000;
    Plan p_always;
    p_always.on = true; p_always.start = 0, p_always.interval = 1;

    sep_ftp
      .set_species({EL,PO,IO})
      .set_probability(0.1)
      .set_plan(p_onetime)
      .set_marker([](auto& p){p.set(flag::_5);});
    ycloud
      .set_species({EL,PO,IO})
      .set_probability(0.1)
      .set_plan(p_onetime)
      .set_marker([](auto& p){p.set(flag::_6);});
    two_small_plasmoids
      .set_species({EL,PO,IO})
      .set_probability(0.1)
      .set_plan(p_onetime)
      .set_marker([](auto& p){p.set(flag::_7);});
    one_big_plasmoid
      .set_species({EL,PO,IO})
      .set_probability(0.1)
      .set_plan(p_onetime)
      .set_marker([](auto& p){p.set(flag::_8);});
    black_region
      .set_species({EL,PO})
      .set_probability(0.1)
      .set_plan(p_onetime)
      .set_marker([](auto& p){p.set(flag::_9);});
    backflow
      .set_species({EL})
      .set_probability(1.01)
      .set_plan(p_always)
      .set_marker
      ([](auto&p) {
         constexpr auto f // encoding of radius
           = [](real_t r) {
               return static_cast<int>( std::min( (r-6.0_r)*0.5_r, 0.001_r ) );
             };
         // use one bit to flag, use three bits to encode farthest distance traveled
         int r = f(std::exp(p.q(0)));
         int r_max = p.is(flag::_11) + p.is(flag::_12)*2 + p.is(flag::_13)*4;
         (r < r_max and r < f(12.01_r)) ? p.set(flag::_10) : p.reset(flag::_10);
         r_max = std::max(r_max,r);
         (r_max & 1) ? p.set(flag::_11) : p.reset(flag::_11);
         (r_max & 2) ? p.set(flag::_12) : p.reset(flag::_12);
         (r_max & 4) ? p.set(flag::_13) : p.reset(flag::_13);
       });
    grand_tot.set_species({EL,PO,IO})
      .set_probability(1.01)
      .set_plan(save_tracing_plan)
      .set_marker
      ([](auto&p) {
         if ( p.is(flag::_5) or p.is(flag::_6) or p.is(flag::_7)
              or p.is(flag::_8) or p.is(flag::_9) or p.is(flag::_10) ) {
           pic::trace(p);
         } else {
           pic::untrace(p);
         }
       });

    pus.emplace_back(no_ph.Clone());
    pus.emplace_back(annih.Clone());
    pus.emplace_back(escape.Clone());
    pus.emplace_back(pu.Clone());
    pus.emplace_back(atm.Clone());
    pus.emplace_back(asym_lo.Clone());
    pus.emplace_back(asym_hi.Clone());
    pus.emplace_back(analyzer.Clone());
    // FIXME migrate need more memory check
    pus.emplace_back(migrate.Clone()); // After this line, particles are all within borders.

    if ( save_tracing_plan.on ) {
      pus.emplace_back(sep_ftp.Clone());
      pus.emplace_back(ycloud.Clone());
      pus.emplace_back(two_small_plasmoids.Clone());
      pus.emplace_back(one_big_plasmoid.Clone());
      pus.emplace_back(black_region.Clone());
      pus.emplace_back(backflow.Clone());
      pus.emplace_back(grand_tot.Clone());
    }

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
          B[1](I) = B_theta_star( grid[0].absc(I[0], 0.5 * B[1].offset()[0]), grid[1].absc(I[1], 0.5 * B[1].offset()[1]), 0, 0 );
        }
      }

    } ic;
    ic[0] = { 0, supergrid[0].dim() };
    ic[1] = { 0, supergrid[1].dim() + 1 }; // NOTE +1 to include upper boundary

    return ic;
  }
}

#include "io/exportee.hpp"
#include "msh/mesh_shape_interplay.hpp"
#include <cassert>

namespace pic {
  constexpr bool is_collinear_mesh = false; // FIXME this is an ad hoc fix

  void export_prior_hook( const map<PtcArray>& particles, const map<Properties>& properties,
                          const Field<3>& E, const Field<3>& B, const JField& J,  const Grid& grid, const Ensemble& ens,
                          real_t dt, int timestep ) {
    { // pair creation counter
      auto& pc = RTD::data().pc_counter;
      ens.reduce_to_chief( mpi::by::SUM, pc[0].data().data(), pc[0].data().size() );
    }

    { // skin depth
      auto& skin_depth = RTD::data().skin_depth;
      skin_depth = {J.mesh()};
      skin_depth.reset();

      for ( auto sp : particles ) {
        auto q2m = properties[sp].charge_x * properties[sp].charge_x / properties[sp].mass_x;
        for ( const auto& ptc : particles[sp] ) {
          if ( !ptc.is(flag::exist) ) continue;
          Index I;
          for ( int i = 0; i < DGrid; ++i )
            I[i] = grid[i].csba( ptc.q(i) );
          skin_depth[0](I) += q2m * ptc.frac();
        }
      }
      ens.reduce_to_chief( mpi::by::SUM, skin_depth[0].data().data(), skin_depth[0].data().size() );
      if ( ens.is_chief() ) {
        for ( const auto& I : apt::Block(apt::range::begin(skin_depth.mesh().range()), apt::range::end(skin_depth.mesh().range())) ) {
          real_t r = grid[0].absc(I[0], 0.5);
          real_t theta = grid[1].absc(I[1], 0.5);
          real_t h = Metric::h<2>(r,theta) / (wpic2 * grid[0].delta() * grid[0].delta());
          auto& v = skin_depth[0](I);
          v = std::sqrt(h / v);
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

  void average_when_downsampled ( IOField& fds, const IOGrid& , int num_comps, const mpi::CartComm& ) {
    int factor = POW(export_plan.downsample_ratio, pic::DGrid);
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

  apt::array<real_t,3> pair_creation_rate ( Index I, const Grid& grid, const Field<3>& ,
                                            const Field<3>& , const JField&  ) {
    auto x = msh::interpolate( RTD::data().pc_counter, I2std(I), ShapeF() );
    return { x[0] / RTD::data().pc_cumulative_time, 0, 0};
  }

  template < particle::species SP >
  apt::array<real_t,3> frac_J_sp ( Index I, const Grid& grid, const Field<3>& ,
                                   const Field<3>& , const JField& J ) {
    auto q = I2std(I);
    auto j_sp = msh::interpolate( RTD::data().Jsp[SP], q, ShapeF() );
    auto j = msh::interpolate( J, q, ShapeF() );
    for ( int i = 0; i < 3; ++i ) {
      if ( j[i] == 0 ) j_sp[i] = 0;
      else j_sp[i] /= j[i];
    }
    return { real_t(j_sp[0]), real_t(j_sp[1]), real_t(j_sp[2]) };
  }

  apt::array<real_t,3> skin_depth ( Index I, const Grid& grid, const Field<3>& ,
                                    const Field<3>& , const JField& ) {
    return { msh::interpolate( RTD::data().skin_depth, I2std(I), ShapeF() )[0], 0, 0  };
  }

  auto set_up_field_export() {
    std::vector<::io::FieldExportee<real_export_t, DGrid, real_t, real_j_t>*> fexps;
    {
      using FA = ::io::FexpTbyFunction<real_export_t, DGrid, real_t, real_j_t>;

      fexps.push_back( new FA ( "E", 3, field_self<0>, average_when_downsampled) );
      fexps.push_back( new FA ( "B", 3, field_self<1>, average_when_downsampled) );
      fexps.push_back( new FA ( "J", 3, field_self<2>, average_and_divide_flux_by_area) );
      fexps.push_back( new FA ( "PairCreationRate", 1, pair_creation_rate, nullptr) );
      fexps.push_back( new FA ( "SkinDepth", 1, skin_depth, average_when_downsampled) );

      if ( RTD::data().is_export_Jsp ) {
        using namespace particle;
        for ( auto sp : RTD::data().Jsp ) {
          switch(sp) {
          case species::electron :
            fexps.push_back( new FA ( "fJ_Electron", 3, frac_J_sp<species::electron>, average_when_downsampled) ); break;
          case species::positron :
            fexps.push_back( new FA ( "fJ_Positron", 3, frac_J_sp<species::positron>, average_when_downsampled) ); break;
          case species::ion :
            fexps.push_back( new FA ( "fJ_Ion", 3, frac_J_sp<species::ion>, average_when_downsampled) ); break;
          default: ;
          }
        }
      }
    }
    return fexps;
  }

}

namespace pic {
  apt::array<real_t,3> ptc_num ( const Properties& prop, const typename PtcArray::const_particle_type& ptc ) {
    return { 1.0, 0.0, 0.0 };
  }

  apt::array<real_t,3> ptc_energy ( const Properties& prop, const typename PtcArray::const_particle_type& ptc ) {
    return { std::sqrt( (prop.mass_x > 0.01) + apt::sqabs(ptc.p()) ), 0.0, 0.0 };
  }

  apt::array<real_t,3> ptc_momentum ( const Properties& prop, const typename PtcArray::const_particle_type& ptc ) {
    return { ptc.p(0), ptc.p(1), ptc.p(2) };
  }

  auto set_up_particle_export() {
    std::vector<::io::PtcExportee<real_export_t, DGrid, real_t, Specs>*> pexps;
    {
      using PA = ::io::PexpTbyFunction<real_export_t,DGrid,real_t,Specs>;
      pexps.push_back( new PA ("Num", 1, ptc_num, nullptr) );

      pexps.push_back( new PA ("E", 1, ptc_energy, nullptr) );

      pexps.push_back( new PA ("P", 3, ptc_momentum, nullptr) );
    }

    return pexps;
  }
}

namespace pic {
  void export_post_hook() {
    auto& pc = RTD::data().pc_counter;
    std::fill( pc[0].data().begin(), pc[0].data().end(), 0 );
    RTD::data().pc_cumulative_time = 0;
    if ( RTD::data().is_export_Jsp ) {
      // clear Jsp to save some space
      for ( auto sp : RTD::data().Jsp )
        RTD::data().Jsp[sp] = {};
    }
    RTD::data().skin_depth = {};
  }
}

#include <sstream>
#include "apt/print.hpp"

namespace pic {
  std::string proofread(std::string indent) {
    std::ostringstream o;
    real_t gamma_0 = Omega * Omega * mu;
    o << indent << "gamma_0=" << apt::fmt("%.0f", gamma_0 ) << std::endl;
    o << indent << "Np=" << apt::fmt("%.1f", 2 * Omega * mu / wpic2 ) << std::endl;
    o << indent << "(w_pic dt)^2 = " << apt::fmt("%.4f", wpic2 * dt * dt ) << std::endl;
    o << indent << "re=" << apt::fmt("%.4f", r_e() ) << std::endl;
    o << indent << "Ndot_GJ=" << apt::fmt("%.4e", gamma_0 / ( 180 * r_e() ) ) << std::endl;

    {
      using namespace particle;
      o << indent << "ATM: N_atm_floor=" << apt::fmt("%.1f", N_atm_floor());
      o << ", atm_x=" << apt::fmt("%.1f", atm_x);
      o << ", v_th=" << apt::fmt("%.2f", v_th)
        << ", g=" << apt::fmt("%.2f", gravity_strength) << std::endl;
    }

    {
      using namespace particle;
      o << indent << "PC: gamma_fd=" << apt::fmt("%.0f", gamma_fd);
      o << ", E_ph=" << apt::fmt("%.0f", E_ph);

      o << indent << "    gamma_RRL=" << gamma_fd * std::pow(gamma_fd / E_ph, 0.5);
      // o << ", L_CR/L_sd=" << E_ph * Ndot_fd * Omega * std::pow( gamma_0 / gamma_fd, 3.0 );
    }
    return o.str();
  }
}

#endif
