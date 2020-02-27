#ifndef _PIC_IMPL_HPP_
#define _PIC_IMPL_HPP_

#include "../examples/pic_prior_impl.hpp"

#include "apt/numeric.hpp"

#include "pic/module_range.hpp"
#include "pic/forces/gravity.hpp"
#include "pic/forces/landau0.hpp"

#include "metric/log_spherical.hpp"

#include "field/old_field_solver/updater.hpp" // FIXME

#include "io/exportee_by_function.hpp"

namespace pic {
  using Metric = metric::LogSpherical<real_t>;

  constexpr double PI = std::acos(-1.0);

  inline constexpr const char* project_name = "Pulsar64";
  inline constexpr const char* datadir_prefix = "../Data/";

  inline constexpr apt::array<int,DGrid> dims = { 1, 1 };
  inline constexpr apt::array<bool,DGrid> periodic = {false,false};
  inline constexpr real_t dt = 0.01;

  constexpr Grid supergrid
  = {{ { 0.0, std::log(30.0), 64 }, { 0.0, PI, 64 } }};

  inline constexpr real_t wdt_pic = 1.0 / 30.0;
  inline constexpr real_t mu = 3750; // set the impact of field strength on particle

  inline constexpr int initial_timestep = 0;
  inline constexpr int total_timesteps = 100;

  constexpr real_t classic_electron_radius () noexcept {
    real_t res = wdt_pic * wdt_pic / ( 4.0 * std::acos(-1.0l) * dt * dt);
    apt::foreach<0,DGrid>( [&res](const auto& g) { res *= g.delta(); }, supergrid );
    return res;
  }
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
  inline constexpr std::optional<int> msperf_max_entries {100};
  inline constexpr auto msperf_qualified =
    []( const std::optional<dye::Ensemble<DGrid>>& ens_opt ) -> bool { return true; };

  inline constexpr ModuleRange vitals_mr { true, 0, 100 };

  inline constexpr int cout_ts_interval = 100;

  inline constexpr ModuleRange tracing_mr { false, 1, 1000 };
  inline constexpr int num_tracing_parts = 4;
}

namespace pic {
  inline constexpr real_t Omega = 1.0 / 6.0;

  constexpr int star_interior = 5;

  constexpr int myguard = std::max(1, ( pic::ShapeF::support() + 3 ) / 2 ); // NOTE minimum number of guards of J on one side is ( supp + 3 ) / 2

  constexpr real_t omega_spinup ( real_t time ) noexcept {
    return std::min<real_t>( time / 4.0, 1.0 ) * Omega;
  }

  constexpr real_t B_r_star( real_t lnr, real_t theta, real_t , real_t time ) noexcept {
    return 2.0 * std::cos(theta) * std::exp(-3.0 * lnr);
  }
  constexpr real_t B_theta_star( real_t lnr, real_t theta, real_t , real_t time ) noexcept {
    return std::sin(theta) * std::exp(-3.0 * lnr);
  }
  constexpr real_t B_phi_star( real_t lnr, real_t theta, real_t , real_t time ) noexcept {
    return 0;
  }

  constexpr real_t E_r_star( real_t lnr, real_t theta, real_t , real_t time ) noexcept {
    return omega_spinup(time) * std::exp(- 2 * lnr) * std::sin( theta ) * std::sin(theta);
  }
  constexpr real_t E_theta_star( real_t lnr, real_t theta, real_t , real_t time ) noexcept {
    return - omega_spinup(time) * std::exp(- 2 * lnr ) * std::sin( 2*theta );
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
  struct RTD {
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
      is_export_Jsp = true;
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
  namespace yee = ::field::yee;
  template < ::field::offset_t Ftype >
  constexpr auto Diff = ::field::Diff<DGrid,real_t,Ftype>;
  constexpr auto HHk = ::field::HHk<real_t>;
  constexpr auto diff_zero = ::field::diff_zero<DGrid,real_t>;
  constexpr auto diff_one = ::field::diff_one<real_t>;

  auto set_up_field_actions() {
    std::vector<std::unique_ptr<FieldAction>> fus;
    namespace range = apt::range;

    ::field::OldSolve<real_t,DGrid,real_j_t> fu;
    {
      const int guard = fu.guard();
      fu.setName("OldSolve");
      fu[0] = { 0, pic::supergrid[0].dim(), guard };
      fu[1] = { 0, pic::supergrid[1].dim(), guard };

      constexpr auto PREJ = 4.0 * std::acos(-1.0l) * pic::classic_electron_radius();
      fu.set_fourpi(PREJ);
      fu.set_mu(pic::mu);
      fu.set_magnetic_pole(2);
      fu.set_damping_rate(10.0);
      fu.set_surface_indent(5);
      fu.set_damp_indent(43);
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
                            v_g = ( &v_g == &v_b ) ? 0.0 : - v_b;
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
  constexpr real_t gravity_strength = 0.5;
  real_t landau0_B_thr = 0.1;

  constexpr real_t gamma_fd = 20.0;
  constexpr real_t gamma_off = 15.0;
  constexpr real_t Ndot_fd = 0.25;
  constexpr real_t E_ph = 4.0;

  auto set_up_particle_properties() {
    map<Properties> properties;
    {
      properties.insert(species::electron, {1,-1,"electron","el"});
      //properties.insert(species::positron, {1,1,"positron","po"});
      properties.insert(species::ion, { 5, 1, "ion","io"});
      //properties.insert(species::photon, { 0, 0, "photon","ph" });
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
        // force.add( landau0, landau0_B_thr );

        force.Register(sp);
      }
      if ( properties.has(species::positron) ) {
        auto sp = species::positron;
        Force force;
        const auto& prop = properties[sp];

        force.add( lorentz, prop.charge_x / prop.mass_x );
        force.add( gravity, gravity_strength );
        // force.add( landau0, landau0_B_thr );

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

      ep_scat.eligs.push_back([](const Ptc_t& ptc){ return ptc.q(0) < std::log(9.0); });

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
        // ep_scat.Register( species::electron );
        // ep_scat.Register( species::positron );
      }
    }

    if ( properties.has(species::photon) ) {
      ::particle::Scat<real_t,Specs> photon_scat;
      // Photons are free to roam across all domain. They may produce pairs outside light cylinder
      photon_scat.eligs.push_back([](const Ptc_t& ptc) { return true; });
      scat::MagneticConvert<real_t,Specs>::B_thr = 0.1;
      scat::MagneticConvert<real_t,Specs>::mfp = 0.2;
      photon_scat.channels.push_back( scat::MagneticConvert<real_t,Specs>::test );

      scat::TwoPhotonCollide<real_t,Specs>::mfp = 5.0;
      photon_scat.channels.push_back( scat::TwoPhotonCollide<real_t,Specs>::test );

      photon_scat.impl = scat::PhotonPairProduction<real_t,Specs>;

      photon_scat.Register( species::photon );
    }

    return properties;
  }

}

namespace pic {
  constexpr real_t N_atm_floor = std::exp(1.0) * 2.0 * Omega * mu * std::pow( dt / wdt_pic, 2.0 );
  constexpr real_t N_atm_x = 1.0; // TODO
  constexpr real_t v_th = 0.2;

  using R = real_t;

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
                        R dt, int timestep, util::Rng<R>& rng
                        ) override {
        if ( !RTD::data().is_export_Jsp || !export_data_mr.is_do(timestep) ) {
          _pu( particles,J, new_ptc_buf, properties, E, B, grid, ens, dt, timestep, rng );
        } else {
          auto& Jsp = RTD::data().Jsp;
          // store J by species separately for data export
          for ( auto sp : particles ) {
            map<PtcArray> ptcs_sp;
            ptcs_sp.insert(sp);
            std::swap( ptcs_sp[sp], particles[sp] ); // FIXME make sure there is no copying
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
      pu0.set_update_q(Metric::geodesic_move<apt::vVec<R,3>, apt::vVec<R,3>>);

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
      R _v_th = 0.0;
      R _N_atm = 0.0;
      R _min_frac = 1e-6; // over fracs larger than this will be injected
      R (*_omega_t) ( R time ) = nullptr;

    public:
      auto& set_thermal_velocity(R v) { _v_th = v; return *this; }
      auto& set_number_in_atmosphere(R N) { _N_atm = N; return *this; }
      auto& set_minimal_fraction( R x ) { _min_frac = x; return *this; }
      auto& set_normal_direction( int n ) { _n = n; return *this; }
      auto& set_omega_t(R (*omega_t) ( R )) { _omega_t = omega_t; return *this; }
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
                        R dt, int timestep, util::Rng<R>& rng
                        ) override {
        _count_n.resize( {apt::make_range(range::begin(*this),range::end(*this),0)} );
        _count_p.resize( {apt::make_range(range::begin(*this),range::end(*this),0)} );

        apt::array<R,DGrid> lb;
        apt::array<R,DGrid> ub;
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
                  idx[i] = ( x.q(i) - grid[i].lower() ) / grid[i].delta(); // NOTE used grid.lower instead of lb, important
                count[0](idx) += x.frac(); // add by fraction
              }
            };

        f_count( _count_n, particles[_negaon] );
        f_count( _count_p, particles[_posion] );

        { // parallelizae TODO optimize
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
                nB[0] = 0.5 * ( B[0][li] + B[0][li + m.stride(1)] );
                nB[1] = 0.5 * ( B[1][li] + B[1][li + m.stride(0)] );
                nB[2] = 0.25 * ( B[2][li] + B[2][li + m.stride(0)] + B[2][li + m.stride(1)] + B[2][li + m.stride(0) + m.stride(1)] );
              } else if (DGrid == 3){
              nB[0] = 0.25 * ( B[0][li] + B[0][li + m.stride(1)] + B[0][li + m.stride(2)] + B[0][li + m.stride(1) + m.stride(2)] );
              nB[1] = 0.25 * ( B[1][li] + B[1][li + m.stride(2)] + B[1][li + m.stride(0)] + B[1][li + m.stride(2) + m.stride(0)] );
              nB[2] = 0.25 * ( B[2][li] + B[2][li + m.stride(0)] + B[2][li + m.stride(1)] + B[2][li + m.stride(0) + m.stride(1)] );
            }
            if ( apt::abs(nB) == 0.0 ) nB = {1.0, 0.0, 0.0}; // use radial direction as default
            else nB /= apt::abs(nB);
          }

          Vec3 p{};
          p[2] = _omega_t( timestep * dt ) * std::exp(q[0]) * std::sin(q[1]); // corotating

          // replenish
          R quota = _N_atm * std::sin(q[1]) - N_pairs;
          while ( quota > _min_frac ) {
            auto q_ptc = q;
            R frac = std::min( (R)1.0, quota );
            quota -= (R)1.0;

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

      atm.set_thermal_velocity(v_th).set_number_in_atmosphere(N_atm_x * N_atm_floor);
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
                                R, int, util::Rng<R>&) override {
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
                                R dt, int timestep, util::Rng<R>& rng) override {
        // Put particles where they belong after scattering
        assert(new_ptc_buf != nullptr);

        bool is_trace_new = ( timestep > 0 and ( timestep % 100 ) == 0 );

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
                I[j] = ( ptc.q(j) - grid[j].lower() ) / grid[j].delta();
              RTD::data().pc_counter[0](I) += ptc.frac();

              // trace electrons near Y point
              if ( is_trace_new
                   and std::log(6.0) < ptc.q(0) and ptc.q(0) < std::log(7.0)
                   and 1.47 < ptc.q(1) and ptc.q(1) < 1.67
                   and rng.uniform() < 0.01 ) {
                Tracer::trace(ptc);
              }
            }
          }

          particles[this_sp].push_back(std::move(ptc));
        }
        new_ptc_buf->resize(0);
        RTD::data().pc_cumulative_time += dt;

        if ( is_trace_new) {
          // trace primary electrons near foot of upper separatrix.
          for ( auto ptc : particles[species::electron] ) { // TODOL semantics
            if ( not ptc.is(flag::exist) ) continue;

            if ( std::log(1.0) < ptc.q(0)
                 and (35.0*PI / 180.0) < ptc.q(1) and ptc.q(1) < (45.0*PI/180.0)
                 and rng.uniform() < 0.01 ) {
              Tracer::trace(ptc);
            }
          }
        }
      }
    } analyzer;
    {
      analyzer.setName("NewPtcAnalysis");
    }

    pus.emplace_back(pu.Clone());
    pus.emplace_back(atm.Clone());
    pus.emplace_back(asym_lo.Clone());
    pus.emplace_back(asym_hi.Clone());
    pus.emplace_back(analyzer.Clone());
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
        auto q2m = properties[sp].charge_x * properties[sp].charge_x / static_cast<R>(properties[sp].mass_x);
        for ( const auto& ptc : particles[sp] ) {
          if ( !ptc.is(flag::exist) ) continue;
          Index I;
          for ( int i = 0; i < DGrid; ++i )
            I[i] = ( ptc.q(i) - grid[i].lower() ) / grid[i].delta();
          skin_depth[0](I) += q2m * ptc.frac();
        }
      }
      ens.reduce_to_chief( mpi::by::SUM, skin_depth[0].data().data(), skin_depth[0].data().size() );
      if ( ens.is_chief() ) {
        for ( const auto& I : apt::Block(apt::range::begin(skin_depth.mesh().range()), apt::range::end(skin_depth.mesh().range())) ) {
          R r = grid[0].absc(I[0], 0.5);
          R theta = grid[1].absc(I[1], 0.5);
          R h = Metric::h<2>(r,theta) * dt * dt / ( wdt_pic * wdt_pic * grid[0].delta() * grid[0].delta() );
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

  void divide_flux_by_area ( IOField& fds, const IOGrid& grid, int num_comps, const mpi::CartComm& ) {
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

      fexps.push_back( new FA ( "E", 3, field_self<0>, nullptr) );
      fexps.push_back( new FA ( "B", 3, field_self<1>, nullptr) );
      fexps.push_back( new FA ( "J", 3, field_self<2>, divide_flux_by_area) );
      fexps.push_back( new FA ( "PairCreationRate", 1, pair_creation_rate, nullptr) );
      fexps.push_back( new FA ( "SkinDepth", 1, skin_depth, nullptr) );

      if ( RTD::data().is_export_Jsp ) {
        using namespace particle;
        for ( auto sp : RTD::data().Jsp ) {
          switch(sp) {
          case species::electron :
            fexps.push_back( new FA ( "fJ_Electron", 3, frac_J_sp<species::electron>, nullptr) ); break;
          case species::positron :
            fexps.push_back( new FA ( "fJ_Positron", 3, frac_J_sp<species::positron>, nullptr) ); break;
          case species::ion :
            fexps.push_back( new FA ( "fJ_Ion", 3, frac_J_sp<species::ion>, nullptr) ); break;
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
    return { std::sqrt( (prop.mass_x != 0) + apt::sqabs(ptc.p()) ), 0.0, 0.0 };
  }

  apt::array<real_t,3> ptc_momentum ( const Properties& prop, const typename PtcArray::const_particle_type& ptc ) {
    return { ptc.p(0), ptc.p(1), ptc.p(2) };
  }

  auto set_up_particle_export() {
    std::vector<::io::PtcExportee<real_export_t, DGrid, real_t, Specs>*> pexps;
    {
      using PA = ::io::PexpTbyFunction<real_export_t,DGrid,real_t,Specs,particle::induced_shapef_t<ShapeF,downsample_ratio>>;
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

  std::string characteristics(std::string indent) {
    std::ostringstream o;
    auto gamma_0 = std::pow(Omega,2.0) * mu;
    o << indent << "gamma_0=" << apt::fmt("%.0f", gamma_0 ) << std::endl;
    o << indent << "w_pic dt=" << apt::fmt("%.4f", wdt_pic ) << std::endl;
    o << indent << "re=" << apt::fmt("%.4f", classic_electron_radius() ) << std::endl;
    o << indent << "Ndot_GJ=" << apt::fmt("%.4e", gamma_0 / ( 2 * PI * classic_electron_radius() ) ) << std::endl;

    {
      using namespace particle;
      o << indent << "ATM: N_atm_floor=" << apt::fmt("%.1f", N_atm_floor);
      o << ", multiplicity=" << apt::fmt("%.1f", N_atm_x);
      o << ", v_th=" << apt::fmt("%.2f", v_th) << ", g=" << apt::fmt("%.2f", gravity_strength) << std::endl;
    }

    {
      using namespace particle;
      o << indent << "PC: gamma_fd=" << apt::fmt("%.0f", gamma_fd);
      o << ", gamma_off=" << apt::fmt("%.0f", gamma_off);
      o << ", E_ph=" << apt::fmt("%.0f", E_ph);
      o << ", Ndot_fd=" << apt::fmt("%.2f", Ndot_fd) << std::endl;

      o << indent << "    gamma_RRL=" << gamma_fd * std::pow(gamma_fd / E_ph, 0.5);
      o << ", L_CR/L_sd=" << E_ph * Ndot_fd * Omega * std::pow( gamma_0 / gamma_fd, 3.0 );
    }
    return o.str();
  }
}

#endif
