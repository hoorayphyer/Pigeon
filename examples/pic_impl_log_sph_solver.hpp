#ifndef _PIC_IMPL_HPP_
#define _PIC_IMPL_HPP_

#include "../examples/pic_prior_impl.hpp"
#include "apt/numeric.hpp"
#include "metric/log_spherical.hpp"
#include "field/log_spherical_solver/updater.hpp"
#include "io/exportee_by_function.hpp"

namespace pic {
  constexpr double PI = std::acos(-1.0);

  inline constexpr apt::array<int,DGrid> dims = { 1, 1 };
  inline constexpr apt::array<bool,DGrid> periodic = {false,false};

  constexpr Grid supergrid
  = {{ { 0.0, std::log(30.0), 256 }, { 0.0, PI, 256 } }};

  real_t Omega;
  real_t dt;
  real_t mu;

  real_t damping_thickness;
  real_t damping_rate;
  real_t spinup_time;

  std::string project_name;
  std::string datadir_prefix;
  int total_timesteps;

  void load_configuration(const ConfFile_t& conf) {
    safe_set(dt, conf["dt"]);
    safe_set(mu, conf["mu"]);
    safe_set(Omega, conf["Omega"]);

    safe_set(damping_thickness, conf["damping"]["thickness"]);
    safe_set(damping_rate, conf["damping"]["rate"]);
    safe_set(spinup_time, conf["spinup_time"]);

    project_name = conf["project_name"].value_or("Unnamed"sv);
    datadir_prefix = conf["datadir_prefix"].value_or("../Data/"sv);
    total_timesteps = conf["total_timesteps"].value_or(100);

    print_timestep_to_stdout_interval = conf["plans"]["print_timestep_to_stdout_interval"].value_or(100);

    export_plan.on = conf["plans"]["export"]["on"].value_or(true);
    export_plan.start = conf["plans"]["export"]["start"].value_or(0);
    safe_set(export_plan.interval, conf["plans"]["export"]["interval"]);
    export_plan.num_files = conf["plans"]["export"]["num_files"].value_or(1);
    export_plan.downsample_ratio = conf["plans"]["export"]["downsample_ratio"].value_or(1);

    profiling_plan.on = conf["plans"]["profiling"]["on"].value_or(false);
    profiling_plan.start = conf["plans"]["profiling"]["start"].value_or(0);
    safe_set(profiling_plan.interval, conf["plans"]["profiling"]["interval"]);
    if ( auto n = conf["plans"]["profiling"]["max_entries"].value<int64_t>() ) {
      profiling_plan.max_entries = {*n};
    }
  }

  int damping_begin_index() {
    return (std::log(std::exp(supergrid[0].upper()) - damping_thickness) -
            supergrid[0].lower()) / supergrid[0].delta();
  }
}

namespace pic {
  constexpr int field_op_inv_precision = 4;
  constexpr int myguard = std::max(
      ::field::LogSphericalSolver<real_t, DGrid, real_j_t>::min_guard(field_op_inv_precision, false),
      (pic::ShapeF::support() + 3) / 2); // NOTE minimum number of guards of J
                                         // on one side is ( supp + 3 ) / 2

  real_t omega_spinup ( real_t time ) noexcept {
    return std::min<real_t>( time / spinup_time, 1.0 ) * Omega;
  }

  real_t B_r_star( real_t lnr, real_t theta, real_t , real_t ) noexcept {
    return pic::mu * 2.0 * std::cos(theta) * std::exp(-3.0 * lnr);
  }
  real_t B_theta_star( real_t lnr, real_t theta, real_t , real_t ) noexcept {
    return pic::mu * std::sin(theta) * std::exp(-3.0 * lnr);
  }
  constexpr real_t B_phi_star( real_t, real_t, real_t , real_t ) noexcept {
    return 0;
  }

  real_t E_r_star( real_t lnr, real_t theta, real_t , real_t time ) noexcept {
    real_t sin = std::sin(theta);
    return pic::mu * omega_spinup(time) * std::exp(- 2 * lnr) * sin * sin;
  }
  real_t E_theta_star( real_t lnr, real_t theta, real_t , real_t time ) noexcept {
    return - pic::mu * omega_spinup(time) * std::exp(- 2 * lnr ) * std::sin( 2*theta );
  }
  constexpr real_t E_phi_star( real_t, real_t, real_t , real_t ) noexcept {
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
  struct RTD {
  public:
    static RTD& data() {
      static RTD r;
      return r;
    }

    map<real_t> N_scat {};

    void init( const map<Properties>& properties, const Grid& localgrid ) {
      for ( auto sp : properties )
        N_scat.insert( sp, 0 );
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

    constexpr int surface_indent = 5;

    ::field::LogSphericalSolver<real_t,DGrid,real_j_t> fu;
    {
      const int guard = myguard;
      fu.setName("LogSphericalSolver");
      fu[0] = { surface_indent, pic::supergrid[0].dim(), guard };
      fu[1] = { 0, pic::supergrid[1].dim()+1, guard };

      fu.set_fourpi( 4.0 * std::acos(-1.0l) );
      fu.set_alpha(1.0);
      fu.set_op_inv_precision(field_op_inv_precision);
      fu.set_surface(supergrid[0].absc(surface_indent));
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
      fu[0] = {0, surface_indent};
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
          static real_t r_damp_b = std::exp( pic::supergrid[0].absc(damping_begin_index(),0) );
          lnr = ( std::exp(lnr) - r_damp_b ) / damping_thickness;
          return 0.5 * lnr * lnr;
        };

      fu.setName("DampingLayer");
      fu[0] = { damping_begin_index(), pic::supergrid[0].dim() };
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

      fu_asym_hi.setName("AxissymmetrizeEBHigher");
      fu_asym_hi[0] = { 0, pic::supergrid[0].dim() };
      fu_asym_hi[1] = { pic::supergrid[1].dim(), pic::supergrid[1].dim() + myguard };
      fu_asym_hi.is_upper_axis(true);
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

  auto set_up_post_resume_actions() {
    struct PostResume : public apt::ActionBase<DGrid> {
      PostResume *Clone() const override {return new PostResume(*this);}

      void
      operator()(const Grid &grid, Field<3> &E, Field<3> &B, JField &J,
                 map<PtcArray> &particles,
                 const std::optional<Ensemble> & ens_opt,
                 int resumed_timestep,
                 std::string this_run_dir) const {}
    };
    return PostResume{};
  }
}

#include "io/exportee.hpp"
#include "msh/mesh_shape_interplay.hpp"
#include <cassert>

namespace pic {
  constexpr bool is_collinear_mesh = false; // FIXME this is an ad hoc fix

  void export_prior_hook(const map<PtcArray> &particles,
                         const map<Properties> &properties, const Field<3> &E,
                         const Field<3> &B, const JField &J, const Grid &grid,
                         const std::optional<mpi::CartComm> &cart_opt,
                         const Ensemble &ens, real_t dt, int timestep) {}
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

  auto set_up_field_export() {
    std::vector<::io::FieldExportee<real_export_t, DGrid, real_t, real_j_t>*> fexps;
    {
      using FA = ::io::FexpTbyFunction<real_export_t, DGrid, real_t, real_j_t>;

      fexps.push_back( new FA ( "E", 3, field_self<0>, average_when_downsampled) );
      fexps.push_back( new FA ( "B", 3, field_self<1>, average_when_downsampled) );
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
