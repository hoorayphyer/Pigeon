#ifndef _PIC_IMPL_HPP_
#define _PIC_IMPL_HPP_

#include "../examples/pic_prior_impl.hpp"

#include "apt/numeric.hpp"

#include "pic/modules.hpp"

#include "metric/log_spherical.hpp"

#include "field/old_field_solver/updater.hpp" // FIXME

#include "io/exportee_by_function.hpp"

namespace pic {
  using Metric = metric::LogSpherical<real_t>;

  constexpr double PI = std::acos(-1.0);

  inline constexpr apt::array<int,DGrid> dims = { 1, 1 };
  inline constexpr apt::array<bool,DGrid> periodic = {false,false};

  constexpr Grid supergrid
  = {{ { 0.0, std::log(30.0), 256 }, { 0.0, PI, 256 } }};

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

  ModuleBase mod_annih {};

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

    print_timestep_to_stdout_interval = conf["mods"]["print_timestep_to_stdout_interval"].value_or(100);

    mod_sort_ptcs.on = conf["mods"]["sort"]["on"].value_or(true);
    mod_sort_ptcs.start = conf["mods"]["sort"]["start"].value_or(0);
    safe_set(mod_sort_ptcs.interval, conf["mods"]["sort"]["interval"]);

    mod_export.on = conf["mods"]["export"]["on"].value_or(true);
    mod_export.start = conf["mods"]["export"]["start"].value_or(0);
    safe_set(mod_export.interval, conf["mods"]["export"]["interval"]);
    mod_export.num_files = conf["mods"]["export"]["num_files"].value_or(1);
    mod_export.downsample_ratio = conf["mods"]["export"]["downsample_ratio"].value_or(1);

    mod_checkpoint.on = conf["mods"]["checkpoint"]["on"].value_or(false);
    mod_checkpoint.start = conf["mods"]["checkpoint"]["start"].value_or(1);
    safe_set(mod_checkpoint.interval, conf["mods"]["checkpoint"]["interval"]);
    mod_checkpoint.num_files = conf["mods"]["checkpoint"]["num_files"].value_or(1);
    mod_checkpoint.max_num_checkpoints = conf["mods"]["checkpoint"]["max_num_checkpoints"].value_or(1);

    mod_load_balance.on = conf["mods"]["load_balance"]["on"].value_or(true);
    mod_load_balance.start = conf["mods"]["load_balance"]["start"].value_or(0);
    safe_set(mod_load_balance.interval, conf["mods"]["load_balance"]["interval"]);
    mod_load_balance.target_load = conf["mods"]["load_balance"]["target_load"].value_or(100000);

    mod_profiling.on = conf["mods"]["profiling"]["on"].value_or(false);
    mod_profiling.start = conf["mods"]["profiling"]["start"].value_or(0);
    safe_set(mod_profiling.interval, conf["mods"]["profiling"]["interval"]);
    if ( auto n = conf["mods"]["profiling"]["max_entries"].value<int64_t>() ) {
      mod_profiling.max_entries = {*n};
    }

    mod_vitals.on = conf["mods"]["vitals"]["on"].value_or(true);
    mod_vitals.start = conf["mods"]["vitals"]["start"].value_or(0);
    safe_set(mod_vitals.interval, conf["mods"]["vitals"]["interval"]);

    mod_tracing.on = conf["mods"]["tracing"]["on"].value_or(false);
    mod_tracing.start = conf["mods"]["tracing"]["start"].value_or(0);
    safe_set(mod_tracing.interval, conf["mods"]["tracing"]["interval"]);
    mod_tracing.num_files = conf["mods"]["tracing"]["num_files"].value_or(1);

    mod_annih.on = conf["mods"]["annihilation"]["on"].value_or(false);
    mod_annih.start = conf["mods"]["annihilation"]["start"].value_or(0);
    safe_set(mod_annih.interval, conf["mods"]["annihilation"]["interval"]);
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
    return std::min<real_t>( time / spinup_time, 1.0 ) * Omega;
  }

  real_t B_r_star( real_t lnr, real_t theta, real_t , real_t time ) noexcept {
    return pic::mu * 2.0 * std::cos(theta) * std::exp(-3.0 * lnr);
  }
  real_t B_theta_star( real_t lnr, real_t theta, real_t , real_t time ) noexcept {
    return pic::mu * std::sin(theta) * std::exp(-3.0 * lnr);
  }
  constexpr real_t B_phi_star( real_t lnr, real_t theta, real_t , real_t time ) noexcept {
    return 0;
  }

  real_t E_r_star( real_t lnr, real_t theta, real_t , real_t time ) noexcept {
    return pic::mu * omega_spinup(time) * std::exp(- 2 * lnr) * std::sin( theta ) * std::sin(theta);
  }
  real_t E_theta_star( real_t lnr, real_t theta, real_t , real_t time ) noexcept {
    return - pic::mu * omega_spinup(time) * std::exp(- 2 * lnr ) * std::sin( 2*theta );
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
  auto set_up_particle_properties() {
    map<Properties> properties;
    return properties;
  }

}

namespace pic {
  auto set_up_particle_actions() {
    std::vector<std::unique_ptr<PtcAction>> pus;
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
                          real_t dt, int timestep ) {}

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
    int factor = POW(mod_export.downsample_ratio, pic::DGrid); // FIXME test mods.export
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

  auto set_up_field_export() {
    std::vector<::io::FieldExportee<real_export_t, DGrid, real_t, real_j_t>*> fexps;
    {
      using FA = ::io::FexpTbyFunction<real_export_t, DGrid, real_t, real_j_t>;

      fexps.push_back( new FA ( "E", 3, field_self<0>, average_when_downsampled) );
      fexps.push_back( new FA ( "B", 3, field_self<1>, average_when_downsampled) );
      fexps.push_back( new FA ( "J", 3, field_self<2>, average_and_divide_flux_by_area) );

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
  auto set_up_particle_export() {
    std::vector<::io::PtcExportee<real_export_t, DGrid, real_t, Specs>*> pexps;
    return pexps;
  }
}

namespace pic {
  void export_post_hook() {}
}

namespace pic {
  std::string proofread(std::string indent) {
    return {};
  }
}

#endif
