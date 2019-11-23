#ifndef _PIC_IMPL_HPP_
#define _PIC_IMPL_HPP_

#include "apt/numeric.hpp"
#include "apt/index.hpp"
#include "apt/grid.hpp"

#include "pic/module_range.hpp"
#include "pic/forces/gravity.hpp"
#include "pic/forces/landau0.hpp"
#include "metric/log_spherical.hpp"

#include "field/field.hpp"

#include "field/old_field_solver/updater.hpp" // FIXME
#include "field/haugbolle_solver/updater.hpp"

#include "particle/properties.hpp"
#include "particle/map.hpp"
#include "particle/array.hpp"
#include "particle/particle.hpp"

#include "dye/ensemble.hpp"

#include "io/exportee_by_function.hpp"

#include "pic.hpp"

namespace pic {
  using Metric = metric::LogSpherical<real_t>;

  constexpr long double PI = std::acos(-1.0l);

  inline constexpr const char* project_name = "Pulsar64";
  inline constexpr const char* datadir_prefix = "../Data/";

  inline constexpr apt::array<int,DGrid> dims = { 1, 1 };
  inline constexpr apt::array<bool,DGrid> periodic = {false,false};
  inline constexpr real_t dt = 0.01;

  constexpr apt::Grid<real_t,DGrid> supergrid
  = {{ { 0.0, std::log(30.0), 64 }, { 0.0, PI, 64 } }};

  inline constexpr real_t wdt_pic = 1.0 / 30.0;
  inline constexpr real_t w_gyro_unitB = 3750; // set the impact of unit field strength on particle

  inline void set_resume_dir( std::optional<std::string>& dir ) {}
  inline constexpr int initial_timestep = 0;
  inline constexpr int total_timesteps = 2000;

  constexpr real_t classic_electron_radius () noexcept {
    real_t res = wdt_pic * wdt_pic / ( 4.0 * std::acos(-1.0l) * dt * dt);
    apt::foreach<0,DGrid>( [&res](const auto& g) { res *= g.delta(); }, supergrid );
    return res;
  }
}

namespace pic {
  inline constexpr ModuleRange sort_particles_mr { true, 0, 100 };

  inline constexpr ModuleRange export_data_mr { true, 0, 50 };
  inline constexpr int pmpio_num_files = 1;
  inline constexpr int downsample_ratio = 1;

  inline constexpr ModuleRange checkpoint_mr { false, 1, 10000 };
  inline constexpr int num_checkpoint_parts = 4;
  inline constexpr int max_num_ckpts = 2;
  inline constexpr std::optional<float> checkpoint_autosave_hourly;

  inline constexpr ModuleRange dlb_mr { false, 1, 1000 };
  inline constexpr std::optional<int (*) ( int )> dlb_init_replica_deploy {}; // take in ensemble label and return the intended number of replicas in that ensemble
  inline constexpr std::size_t dlb_target_load = 100000;

  inline constexpr ModuleRange msperf_mr { true, 0, 10 };
  inline constexpr std::optional<int> msperf_max_entries {};
  inline constexpr auto msperf_qualified =
    []( const std::optional<dye::Ensemble<DGrid>>& ens_opt ) -> bool { return true; };

  inline constexpr ModuleRange vitals_mr { true, 0, 100 };

  inline constexpr int cout_ts_interval = 100;

  inline constexpr ModuleRange tracing_mr { false, 1, 1000 };
  inline constexpr int num_tracing_parts = 4;
}

namespace field {
  template < int DGrid, typename T >
  void axissymmetrize( field::Component<T,DGrid,false> comp, // TODOL semantics on comp
                       void (*f)( decltype(comp[0])& val_guard, decltype(comp[0])& val_bulk ),
                       apt::Index<DGrid> Ib, apt::Index<DGrid> Ie, bool is_upper ) {
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

namespace field {
  using pic::real_t;
  inline constexpr real_t Omega = 1.0 / 6.0;

  constexpr int star_interior = 5;

  constexpr int order_precision = 2; // order precision of update scheme, this means error is O(x^(order + 1));
  constexpr int number_iteration = 4; // number of iterations in inverting the operator in Haugbolle
  static_assert( order_precision % 2 == 0 );
  constexpr auto PREJ = 4.0 * std::acos(-1.0l) * pic::classic_electron_radius() / pic::w_gyro_unitB;

  // FIXME NOTE the 1 + 
  constexpr int myguard = std::max(1 + order_precision * (1+number_iteration) / 2, ( pic::ShapeF::support() + 3 ) / 2 ); // NOTE minimum number of guards of J on one side is ( supp + 3 ) / 2

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

  real_t (*damping_profile)( real_t q ) = nullptr;

  // TODOL annihilation will affect deposition // NOTE one can deposit in the end
}

namespace field {
  constexpr bool use_new = true;

  template < int DGrid, typename R, typename RJ >
  auto set_up_old_field_actions() {
    std::vector<std::unique_ptr<Action<R,DGrid,RJ>>> fus;
    namespace range = apt::range;

    OldSolve<R,DGrid,RJ> fu;
    {
      const int guard = fu.guard();
      fu.setName("OldSolve");
      fu[0] = { 0, pic::supergrid[0].dim(), guard };
      fu[1] = { 0, pic::supergrid[1].dim(), guard };

      fu.set_fourpi(PREJ);
      fu.set_magnetic_pole(2);
      fu.set_damping_rate(10.0);
      fu.set_surface_indent(5);
      fu.set_damp_indent(43);
      fu.set_omega_t(field::omega_spinup);
    }

    struct Axissymmetric : public Action<R,DGrid,RJ> {
    private:
      bool _is_upper_axis = false;
    public:
      auto& is_upper_axis( bool x ) { _is_upper_axis = x; return *this; }
      virtual Axissymmetric* Clone() const { return new Axissymmetric(*this); }
      virtual void operator() ( Field<R,3,DGrid>& E, Field<R,3,DGrid>& B,
                                const Field<RJ,3,DGrid>& , const apt::Grid<R,DGrid>& grid,
                                const mpi::CartComm&, int , R) const override {
        // NOTE Guard cells values are needed when interpolating E and B
        // E_theta, B_r, B_phi are on the axis. All but B_r should be set to zero
        auto assign = []( R& v_g, R& v_b ) noexcept { v_g = v_b; };
        auto neg_assign = []( R& v_g, R& v_b ) noexcept {
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

  template < int DGrid, typename R, typename RJ >
  auto set_up_field_actions() {
    if constexpr ( !use_new ) return set_up_old_field_actions<DGrid,R,RJ>();

    std::vector<std::unique_ptr<Action<R,DGrid,RJ>>> fus;
    namespace range = apt::range;

    Haugbolle<R,DGrid,RJ> fu_bulk;
    {
      auto& fu = fu_bulk;
      fu.setName("Bulk");
      fu[0] = { star_interior + myguard, pic::supergrid[0].dim() - myguard, myguard };
      fu[1] = { myguard, pic::supergrid[1].dim() + 1 - myguard, myguard }; // NOTE +1 because need values on the upper boundary
      fu.set_fourpi(PREJ);
      fu.set_number_iteration(number_iteration);
      fu.set_implicit(0.8); // 0.8 is stable for E_phi FIXME not really, grows after t = 15

      for ( int i = 0; i < 3; ++i ) { // i is coordinate
        for ( int j = 0; j < 3; ++j ) { // j is Fcomp
          if ( i == j ) continue;
          fu.set_D(yee::Etype,j,i, Diff<DGrid,R,yee::Etype>(j,i) );
          fu.set_D(yee::Btype,j,i, Diff<DGrid,R,yee::Btype>(j,i) );
        }
      }
      for ( int Ftype = Ftype; Ftype < 2; ++Ftype )
        for ( int i = 0; i < 3; ++i )
          fu.set_hh( Ftype,i, HHk<R>(i) );
    }

    HaugbolleBdry<R,DGrid,RJ> fu_axis_lo;
    {
      auto& fu = fu_axis_lo;
      fu.setName("LowerAxis");
      fu[0] = fu_bulk[0];
      fu[1] = {0, myguard, {0,myguard}};
      fu.set_boundary(0,-1);

      fu.init_from(fu_bulk);

      apt::Index<DGrid> I = range::begin(fu);
      fu.set_D( yee::Etype, 0, 1, diff_zero<DGrid,R>, I ); // zero because E_r is polar vector, it's derivative sits on the axis. So zero.
      fu.set_D( yee::Etype, 2, 1, diff_axis_Ephi_theta<DGrid,R,false>, I ); // curl E_phi has 1/sin0
      fu.set_D( yee::Etype, 1, 0, diff_zero<DGrid,R>, I ); // not necessary but good to have
      fu.set_D( yee::Btype, 2, 0, diff_zero<DGrid,R>, I );

      fu.set_hh( yee::Etype, 0, 1, diff_one<DGrid,R>, I ); // because numerator will be zero
      fu.set_hh( yee::Etype, 2, 1, diff_one<DGrid,R>, I );
      fu.set_hh( yee::Btype, 2, 0, diff_one<DGrid,R>, I );
    }

    HaugbolleBdry<R,DGrid,RJ> fu_axis_hi;
    {
      auto& fu = fu_axis_hi;
      fu.setName("HigherAxis");
      fu[0] = fu_bulk[0];
      fu[1] = { pic::supergrid[1].dim() + 1 - myguard, pic::supergrid[1].dim() + 1, {myguard,0} }; // NOTE +1 to include upper boundary
      fu.set_boundary(0,1);

      fu.init_from(fu_bulk);

      // ---Axis cells----
      // because we always calculate the derivative that is in the same cell
      // of the field, we really need to look at the high axis to the LOWER
      // bound of a cell that is beyond the axis. The values of derivatives in
      // that cell don't matter unless it's right on the axis. So we first use
      // trivial derivative for all, then only d E_phi / d theta nees a diff_from_left
      apt::Index<DGrid> I = { fu[0].begin(), fu[1].end() -1 };
      for ( int Ftype = 0; Ftype < 2; ++Ftype ) {
        for ( int c = 0; c < 3; ++c ) {
          fu.set_hh( Ftype, c, diff_one<DGrid,R>, I );
          for ( int i = 0; i < DGrid; ++i ) { // NOTE DGrid is just fine
            if ( c == i ) continue; // only need transverse
            fu.set_D( Ftype, c, i, diff_zero<DGrid,R>, I );
          }
        }
      }
      fu.set_D( yee::Etype, 2, 1, diff_axis_Ephi_theta<DGrid,R,true>, I );
    }

    HaugbolleBdry<R,DGrid,RJ> fu_surf;
    {
      auto& fu = fu_surf;
      fu.setName("ConductingSurface");
      fu[0] = { star_interior, star_interior+myguard, {1, myguard} }; // NOTE use 1 because need values from interior
      fu[1] = fu_bulk[1];
      fu.set_boundary( -1, 0 );
      fu.enforce_continuous_transverse_E(true, false);

      fu.init_from(fu_bulk);

      apt::Index<DGrid> I = range::begin(fu);
      fu.set_D( yee::Etype, 1, 0, diff_1sided_from_right<DGrid,R,1,0>, I );
      fu.set_D( yee::Etype, 2, 0, diff_1sided_from_right<DGrid,R,2,0>, I );
    }

    auto set_double_bdry =
      [&fu_bulk]( HaugbolleBdry<R,DGrid,RJ>& fu,
          const HaugbolleBdry<R,DGrid,RJ>& fu_surf,
          const HaugbolleBdry<R,DGrid,RJ>& fu_axis ) {

        fu.setName( fu_surf.name() + fu_axis.name() );
        fu[0] = fu_surf[0];
        fu[1] = fu_axis[1];

        fu.set_boundary( fu_surf.boundary()[0], fu_axis.boundary()[1] );
        fu.enforce_continuous_transverse_E(fu_surf.cont_trans_E()[0], fu_axis.cont_trans_E()[1]);
        fu.init_from(fu_bulk);

        for ( int Ftype = 0; Ftype < 2; ++Ftype ) {
          // r derivative uses fu_surf, theta derivative uses fu_axis
          for ( int C = 1; C < 3; ++C ) {
            for ( const auto& I : apt::Block(range::begin(fu), range::end(fu)) ) {
              fu.set_D( Ftype, (C+0)%3, 0, fu_surf.D(Ftype, (C+0)%3, 0, I), I );
              fu.set_D( Ftype, (C+1)%3, 1, fu_axis.D(Ftype, (C+1)%3, 1, I), I );
            }
          }
          // hh uses that of fu_axis
          for ( int k = 0; k < 3; ++k ) {
            for ( const auto& I : apt::Block(range::begin(fu), range::end(fu)) ) {
              fu.set_hh( Ftype, k, fu_axis.hh(Ftype,k,I), I);
            }
          }
        }
      };


    HaugbolleBdry<R,DGrid,RJ> fu_surf_axis_lo;
    {
      auto& fu = fu_surf_axis_lo;
      set_double_bdry(fu, fu_surf, fu_axis_lo);
      // for r derivative on the axes, just set them to zero
      for ( int i = range::begin(fu, 0); i < range::end(fu,0); ++i ) {
        fu.set_D( yee::Etype, 1, 0, diff_zero<DGrid,R>, {i,fu[1].begin()} );
        fu.set_D( yee::Btype, 2, 0, diff_zero<DGrid,R>, {i,fu[1].begin()} );
      }
    }

    HaugbolleBdry<R,DGrid,RJ> fu_surf_axis_hi;
    {
      auto& fu = fu_surf_axis_hi;
      set_double_bdry(fu, fu_surf, fu_axis_hi);
      // for r derivative on the axes and beyond set them to appropriate diffs so that no values beyond axis are used.
      {
        static_assert(order_precision == 2);
        int last_theta = fu[1].end() - 1;
        for ( int i = range::begin(fu, 0); i < range::end(fu,0); ++i ) {
          fu.set_D(yee::Etype, 1, 0, diff_zero<DGrid,R>, {i,last_theta} );
          fu.set_D(yee::Etype, 2, 0, diff_zero<DGrid,R>, {i,last_theta} ); // beyond axis
          fu.set_D(yee::Btype, 2, 0, diff_zero<DGrid,R>, {i,last_theta} );
          fu.set_D(yee::Btype, 1, 0, diff_zero<DGrid,R>, {i,last_theta} ); // beyond axis
        }
      }

    }

    // NOTE only implemented for LOWER boundary
    // NOTE at lower boundary the surface is at the "+" : |---*+--|---*--->
    struct RotatingConductor : public Action<R,DGrid,RJ> {
    private:
      apt::array<R(*)(R,R,R,R t),3> _E_cond;
      apt::array<R(*)(R,R,R,R t),3> _B_cond;

    public:
      RotatingConductor& set_E( apt::array<R(*)(R,R,R,R t),3> E ) {
        _E_cond = E; return *this;
      }
      RotatingConductor& set_B( apt::array<R(*)(R,R,R,R t),3> B ) {
        _B_cond = B; return *this;
      }
      virtual RotatingConductor* Clone() const { return new RotatingConductor(*this); }
      virtual void operator() ( Field<R,3,DGrid>& E,
                                Field<R,3,DGrid>& B,
                                const Field<RJ,3,DGrid>& ,
                                const apt::Grid<R,DGrid>& grid,
                                const mpi::CartComm&,
                                int timestep,
                                R dt
                                ) const override {
        auto impose =
          [&,this,time=timestep*dt]( auto& F, const auto& F_cond, int comp ) {
            const auto& ofs = F[comp].offset();

            for ( const auto& I : apt::Block(range::begin(*this),range::end(*this)) ) {
              apt::array<R,3> q{};
              for ( int i = 0; i < DGrid; ++i )
                q[i] = grid[i].absc(I[i], ofs[i] * 0.5);
              F[comp](I) = F_cond[comp](q[0],q[1],q[2],time);
            }
          };
        for ( int C = 0; C < 3; ++C )
          impose(E,_E_cond,C);
        for ( int C = 0; C < 3; ++C )
          impose(B,_B_cond,C);
      }
    } fu_cond;
    { // interior
      auto& fu = fu_cond;
      fu.setName("ConductorInterior");
      fu[0] = {0, star_interior};
      fu[1] = {0, pic::supergrid[1].dim()+1};
      fu.require_original_EB(false);
      fu.set_E({ E_r_star, E_theta_star, E_phi_star });
      fu.set_B({ B_r_star, B_theta_star, B_phi_star });
    }

    // NOTE only implemented for UPPER boundary
    struct DampingLayer : public Action<R,DGrid,RJ> {
    private:
      apt::array<R(*)(R,R,R,R t),3> _E_bg;
      apt::array<R(*)(R,R,R,R t),3> _B_bg;
      int _ndir = 0; // normal direction
      R _rate {};
      R(* _profile)(R q_normal) = nullptr;

    public:
      auto& set_E_background( apt::array<R(*)(R,R,R,R t),3> E ) {_E_bg = E; return *this;}
      auto& set_B_background( apt::array<R(*)(R,R,R,R t),3> B ) {_B_bg = B; return *this;}
      auto& set_normal_direction( int n ) {_ndir = n; return *this;}
      auto& set_damping_rate( R rate ) {_rate = rate; return *this;}
      auto& set_damping_profile( R(*p)(R) ) {_profile = p; return *this;}

      virtual DampingLayer* Clone() const { return new DampingLayer(*this); }

      virtual void operator() ( Field<R,3,DGrid>& E,
                                Field<R,3,DGrid>& B,
                                const Field<RJ,3,DGrid>&,
                                const apt::Grid<R,DGrid>& grid,
                                const mpi::CartComm&,
                                int timestep,
                                R dt
                                ) const override {
        const auto transBlock = apt::project_out(_ndir, range::begin(*this), range::end(*this) );

        auto impose =
          [&,this]( auto& F, const auto& F_bg, int comp ) {
            const auto& ofs = F[comp].offset();

            apt::array<R,3> q{};
            for ( int n = range::begin(*this,_ndir); n < range::end(*this,_ndir); ++n ) {
              q[_ndir] = grid[_ndir].absc(n, ofs[_ndir] * 0.5);
              R lambda = 1.0 - _rate * dt * _profile( q[_ndir] );
              for ( const auto& I : transBlock ) {
                for ( int i = 0; i < DGrid; ++i )
                  if ( _ndir != i ) q[i] = grid[i].absc(I[i], ofs[i] * 0.5);

                R f_bg = F_bg[comp](q[0],q[1],q[2],0.0); // time = 0.0, so damp to the initial condition
                F[comp](I) = ( F[comp](I) - f_bg ) * lambda + f_bg;
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
      constexpr R thickness = 7.0; // in normal spherical radius
      constexpr int nb = ( std::log(std::exp(pic::supergrid[0].upper()) - thickness) - pic::supergrid[0].lower() ) / pic::supergrid[0].delta();
      damping_profile =
        [](R lnr) -> R {
          static R r_damp_b = std::exp( pic::supergrid[0].absc(nb,0) );
          lnr = ( std::exp(lnr) - r_damp_b ) / thickness;
          return 0.5 * lnr * lnr;
        };

      fu.setName("DampingLayer");
      fu[0] = { nb, pic::supergrid[0].dim() };
      fu[1] = { 0, pic::supergrid[1].dim() + 1 };
      fu.require_original_EB(false);

      fu.set_normal_direction(0);
      fu.set_damping_rate(10.0);
      fu.set_damping_profile(damping_profile); // TODO check if this works
      fu.set_E_background({ E_r_star, E_theta_star, E_phi_star });
      fu.set_B_background({ B_r_star, B_theta_star, B_phi_star });
    }

    struct Axissymmetric : public Action<R,DGrid,RJ> {
    private:
      bool _is_upper_axis = false;
    public:
      auto& is_upper_axis( bool x ) { _is_upper_axis = x; return *this; }
      virtual Axissymmetric* Clone() const { return new Axissymmetric(*this); }
      virtual void operator() ( Field<R,3,DGrid>& E, Field<R,3,DGrid>& B,
                                const Field<RJ,3,DGrid>& , const apt::Grid<R,DGrid>& grid,
                                const mpi::CartComm&, int , R) const override {
        // NOTE Guard cells values are needed when interpolating E and B
        // E_theta, B_r, B_phi are on the axis. All but B_r should be set to zero
        auto assign = []( R& v_g, R& v_b ) noexcept { v_g = v_b; };
        auto neg_assign = []( R& v_g, R& v_b ) noexcept {
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
      fu_asym_hi[1] = { pic::supergrid[1].dim(), pic::supergrid[1].dim() + myguard - 1 }; // FIXME NOTE -1 to avoid index out of bound on E or B
      fu_asym_hi.is_upper_axis(true);
      fu_asym_hi.require_original_EB(false);
    }

    fus.emplace_back(fu_bulk.Clone());
    fus.emplace_back(fu_surf.Clone());
    fus.emplace_back(fu_axis_lo.Clone());
    fus.emplace_back(fu_axis_hi.Clone());
    fus.emplace_back(fu_surf_axis_lo.Clone());
    fus.emplace_back(fu_surf_axis_hi.Clone());

    fus.emplace_back(fu_cond.Clone());
    fus.emplace_back(fu_damp.Clone());
    fus.emplace_back(fu_asym_lo.Clone());
    fus.emplace_back(fu_asym_hi.Clone());

    return fus;
  }
}

namespace pic {
  template < int DGrid, typename R, typename RJ, template < typename > class S >
  auto set_up_initial_conditions() {
    // local class in a function
    struct InitialCondition : public apt::ActionBase<DGrid> {
      InitialCondition* Clone() const override { return new InitialCondition(*this); }

      void operator() ( const apt::Grid<R,DGrid>& grid,
                        field::Field<R, 3, DGrid>& ,
                        field::Field<R, 3, DGrid>& B,
                        field::Field<RJ, 3, DGrid>& ,
                        particle::map<particle::array<R,S>>&
                        ) const {
        for ( const auto& I : apt::Block(apt::range::begin(*this),apt::range::end(*this)) ) {
          B[0](I) = field::B_r_star( grid[0].absc(I[0], 0.5 * B[0].offset()[0]), grid[1].absc(I[1], 0.5 * B[0].offset()[1]), 0, 0 );
          B[1](I) = field::B_theta_star( grid[0].absc(I[0], 0.5 * B[1].offset()[0]), grid[1].absc(I[1], 0.5 * B[1].offset()[1]), 0, 0 );
        }
      }

    } ic;
    ic[0] = { 0, supergrid[0].dim() };
    ic[1] = { 0, supergrid[1].dim() + 1 }; // NOTE +1 to include upper boundary

    return ic;
  }
}

// FIXME singleton?
namespace {
  template < typename R, int D >
  struct RTD { // runtime data
    static particle::map<R> N_scat;
    static field::Field<R,1,D> pc_counter;
    static R pc_cumulative_time;

    static particle::map<unsigned int> trace_counter;

    static particle::map<field::Field<pic::real_j_t,3,D>> Jsp; // current by species
    static bool is_export_Jsp;

    static field::Field<R,1,D> skin_depth;

    static void init( const particle::map<particle::Properties>& properties, const apt::Grid< R, D >& localgrid ) {
      is_export_Jsp = false;
      for ( auto sp : properties )
        N_scat.insert( sp, 0 );
      if ( is_export_Jsp ) {
        for ( auto sp : properties )
          Jsp.insert( sp, {} );
      }

      apt::Index<D> bulk_dims;
      for ( int i = 0; i < D; ++i ) bulk_dims[i] = localgrid[i].dim();
      auto range = apt::make_range({}, bulk_dims, 0);
      pc_counter = {range};
    };
  };

  template < typename R, int D >
  particle::map<R> RTD<R,D>::N_scat {};

  template < typename R, int D >
  field::Field<R,1,D> RTD<R,D>::pc_counter {};

  template < typename R, int D >
  R RTD<R,D>::pc_cumulative_time = 0;

  template < typename R, int D >
  particle::map<unsigned int> RTD<R,D>::trace_counter {};

  template < typename R, int D >
  particle::map<field::Field<pic::real_j_t,3,D>> RTD<R,D>::Jsp {};

  template < typename R, int D >
  bool RTD<R,D>::is_export_Jsp = false;

  template < typename R, int D >
  field::Field<R,1,D> RTD<R,D>::skin_depth {};
}

namespace particle {
  constexpr pic::real_t N_atm_floor = std::exp(1.0) * 2.0 * field::Omega * pic::w_gyro_unitB * std::pow( pic::dt / pic::wdt_pic, 2.0 );
  constexpr pic::real_t N_atm_x = 1.0; // TODO
  constexpr pic::real_t v_th = 0.2;
  constexpr pic::real_t gravity_strength = 0.5;

  template < int DGrid, typename R, template < typename > class S,
             typename ShapeF, typename RJ >
  auto set_up_particle_actions() {
    namespace range = apt::range;
    std::vector<std::unique_ptr<Action<DGrid,R,S,RJ>>> pus;

    struct Updater_with_J_species: Action<DGrid,R,S,RJ> {
    private:
      Updater<DGrid,R,S,ShapeF,RJ> _pu;

    public:
      void set_updater( Updater<DGrid,R,S,ShapeF,RJ>&& pu) noexcept { _pu = std::move(pu); }

      Updater_with_J_species* Clone() const override { return new Updater_with_J_species(*this); }

      void operator() ( map<array<R,S>>& particles,
                        field::Field<RJ,3,DGrid>& J,
                        std::vector<Particle<R,S>>* new_ptc_buf,
                        const map<Properties>& properties,
                        const field::Field<R,3,DGrid>& E,
                        const field::Field<R,3,DGrid>& B,
                        const apt::Grid< R, DGrid >& grid,
                        const dye::Ensemble<DGrid>* ens,
                        R dt, int timestep, util::Rng<R>& rng
                        ) override {
        if ( !RTD<R,DGrid>::is_export_Jsp || !pic::export_data_mr.is_do(timestep) ) {
          _pu( particles,J, new_ptc_buf, properties, E, B, grid, ens, dt, timestep, rng );
        } else {
          auto& Jsp = RTD<R,DGrid>::Jsp;
          // store J by species separately for data export
          for ( auto sp : particles ) {
            map<array<R,S>> ptcs_sp;
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
      Updater<DGrid,R,S,ShapeF,RJ> pu0;
      pu0.set_update_q(pic::Metric::geodesic_move<apt::vVec<R,3>, apt::vVec<R,3>>);

      pu.setName("MainUpdate");
      pu.set_updater(std::move(pu0));
    }

    Migrator<DGrid,R,S,ShapeF,RJ> migrate;
    {
      migrate.setName("MigrateParticles");
      migrate.set_supergrid(pic::supergrid);
    }

    struct Atmosphere: Action<DGrid,R,S,RJ> {
    private:
      field::Field<R,1,DGrid> _count_n;
      field::Field<R,1,DGrid> _count_p;
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

      void operator() ( map<array<R,S>>& particles,
                        field::Field<RJ,3,DGrid>& J,
                        std::vector<Particle<R,S>>*,
                        const map<Properties>& properties,
                        const field::Field<R,3,DGrid>& E,
                        const field::Field<R,3,DGrid>& B,
                        const apt::Grid< R, DGrid >& grid,
                        const dye::Ensemble<DGrid>* ens,
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
                apt::Index<DGrid> idx;
                for ( int i = 0; i < DGrid; ++i )
                  idx[i] = ( x.q()[i] - lb[i] ) / grid[i].delta();
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
          apt::Vec<R, S<R>::Dim> q{};
          for ( int i = 0; i < DGrid; ++i )
            q[i] = grid[i].absc(I[i], 0.5);

          apt::Vec<R,3> nB;
          { // make nB centered in the cell
            const auto& m = B.mesh();
            auto li = m.linear_index(I);
            if constexpr (DGrid == 2) {
                nB[0] = 0.5 * ( B[0][li] + B[0][li + m.stride()[1]] );
                nB[1] = 0.5 * ( B[1][li] + B[1][li + m.stride()[0]] );
                nB[2] = 0.25 * ( B[2][li] + B[2][li + m.stride()[0]] + B[2][li + m.stride()[1]] + B[2][li + m.stride()[0] + m.stride()[1]] );
              } else if (DGrid == 3){
              nB[0] = 0.25 * ( B[0][li] + B[0][li + m.stride()[1]] + B[0][li + m.stride()[2]] + B[0][li + m.stride()[1] + m.stride()[2]] );
              nB[1] = 0.25 * ( B[1][li] + B[1][li + m.stride()[2]] + B[1][li + m.stride()[0]] + B[1][li + m.stride()[2] + m.stride()[0]] );
              nB[2] = 0.25 * ( B[2][li] + B[2][li + m.stride()[0]] + B[2][li + m.stride()[1]] + B[2][li + m.stride()[0] + m.stride()[1]] );
            }
            if ( apt::abs(nB) == 0.0 ) nB = {1.0, 0.0, 0.0}; // use radial direction as default
            else nB /= apt::abs(nB);
          }

          apt::Vec<R, S<R>::Dim> p{};
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
            *(itr_ne++) = Particle<R,S>( q_ptc, p_ptc, frac, _negaon, birthplace(ens->label()) );
            *(itr_po++) = Particle<R,S>( std::move(q_ptc), std::move(p_ptc), frac, _posion, birthplace(ens->label()) );
          }
        }
      }
    } atm;
    {
      atm.setName("Atmosphere");
      atm[0] = { field::star_interior, field::star_interior + 1 };
      atm[1] = { 0, pic::supergrid[1].dim() };

      atm.set_thermal_velocity(v_th).set_number_in_atmosphere(N_atm_x * N_atm_floor);
      atm.set_positive_charge(species::ion).set_negative_charge(species::electron);
      atm.set_omega_t(field::omega_spinup).set_normal_direction(0);
    }

    struct Axissymmetric : public Action<DGrid,R,S,RJ> {
    private:
      bool _is_upper_axis = false;
    public:
      auto& is_upper_axis( bool x ) { _is_upper_axis = x; return *this; }
      virtual Axissymmetric* Clone() const { return new Axissymmetric(*this); }

      virtual void operator() ( map<array<R,S>>&,
                                field::Field<RJ,3,DGrid>& J,
                                std::vector<Particle<R,S>>*,
                                const map<Properties>&,
                                const field::Field<R,3,DGrid>&,
                                const field::Field<R,3,DGrid>&,
                                const apt::Grid< R, DGrid >& grid,
                                const dye::Ensemble<DGrid>* ,
                                R, int, util::Rng<R>&
                                ) override {
        auto add_assign =
          []( RJ& a, RJ& b ) noexcept {
            a += b;
            b = a;
          };

        auto sub_assign =
          []( RJ& a, RJ& b ) noexcept {
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
      asym_lo[0] = { 0, pic::supergrid[0].dim() };
      asym_lo[1] = { -field::myguard, 1 }; // NOTE +1 so as to set values right on axis
      asym_lo.is_upper_axis(false);

      asym_hi.setName("AxissymmetrizeJHigher");
      asym_hi[0] = { 0 , pic::supergrid[0].dim() };
      asym_hi[1] = { pic::supergrid[1].dim(), pic::supergrid[1].dim() + field::myguard };
      asym_hi.is_upper_axis(true);
    }

    struct ScatteringAnalyzer : public Action<DGrid,R,S,RJ> {
      virtual ScatteringAnalyzer* Clone() const { return new ScatteringAnalyzer(*this); }

      virtual void operator() ( map<array<R,S>>& particles,
                                field::Field<RJ,3,DGrid>& J,
                                std::vector<Particle<R,S>>* new_ptc_buf,
                                const map<Properties>&,
                                const field::Field<R,3,DGrid>&,
                                const field::Field<R,3,DGrid>&,
                                const apt::Grid< R, DGrid >& grid,
                                const dye::Ensemble<DGrid>* ,
                                R dt, int, util::Rng<R>& rng
                                ) override {
        auto& buf = *new_ptc_buf;
        // Put particles where they belong after scattering
        for ( int i = 0; i < buf.size(); ++i ) {
          auto this_sp = buf[i].template get<species>();

          // log scattering events
          RTD<R,DGrid>::N_scat[this_sp] += buf[i].frac();

          // log pair creation locations
          if ( species::electron == this_sp ) {
            apt::Index<DGrid> I; // domain index, not the global index
            for ( int j = 0; j < DGrid; ++j )
              I[j] = ( buf[i].q()[j] - grid[j].lower() ) / grid[j].delta();
            RTD<R,DGrid>::pc_counter[0](I) += 1;

            // trace electrons near Y point
            if ( std::log(6.0) < buf[i].q()[0] and buf[i].q()[0] < std::log(7.0)
                 and 1.47 < buf[i].q()[1] and buf[i].q()[1] < 1.67
                 and rng.uniform() < 0.01 ) {
              buf[i].set(flag::traced);
              buf[i].set(serial_number(RTD<R,DGrid>::trace_counter[species::electron]++));
            }
          }

          particles[this_sp].push_back(std::move(buf[i]));
        }
        buf.resize(0);
        RTD<R,DGrid>::pc_cumulative_time += dt;
      }
    } scat_anlz;
    {
      scat_anlz.setName("ScatteringAnalysis");
    }

    struct ExportPrep : public Action<DGrid,R,S,RJ> {
    public:
      virtual ExportPrep* Clone() const { return new ExportPrep(*this); }

      virtual void operator() ( map<array<R,S>>& particles,
                                field::Field<RJ,3,DGrid>& J,
                                std::vector<Particle<R,S>>*,
                                const map<Properties>& properties,
                                const field::Field<R,3,DGrid>& E,
                                const field::Field<R,3,DGrid>& B,
                                const apt::Grid< R, DGrid >& grid,
                                const dye::Ensemble<DGrid>* ens,
                                R dt, int timestep, util::Rng<R>&
                                ) override {
        if ( !pic::export_data_mr.is_do(timestep) ) return;
        auto& skin_depth = RTD<R,DGrid>::skin_depth;
        skin_depth = {J.mesh()};
        skin_depth.reset();

        for ( auto sp : particles ) {
          auto q2m = std::pow<R>( properties[sp].charge_x, 2.0 ) / R(properties[sp].mass_x);
          for ( const auto& ptc : particles[sp] ) {
            if ( !ptc.is(flag::exist) ) continue;
            apt::Index<DGrid> I;
            for ( int i = 0; i < DGrid; ++i )
              I[i] = ( ptc.q()[i] - grid[i].lower() ) / grid[i].delta();
            skin_depth[0](I) += q2m * ptc.frac();
          }
        }
        ens->reduce_to_chief( mpi::by::SUM, skin_depth[0].data().data(), skin_depth[0].data().size() );
        if ( ens->is_chief() ) {
          for ( const auto& I : apt::Block(apt::range::begin(skin_depth.mesh().range()), apt::range::end(skin_depth.mesh().range())) ) {
            R r = grid[0].absc(I[0], 0.5);
            R theta = grid[1].absc(I[1], 0.5);
            R h = std::sqrt( pic::Metric::h<2>(r,theta) ) * dt / ( pic::wdt_pic * grid[0].delta() );
            auto& v = skin_depth[0](I);
            v = h / v;
          }
        }

      }
    } export_prep;
    {
      export_prep.setName("ExportPrep");
    }

    pus.emplace_back(pu.Clone());
    pus.emplace_back(atm.Clone());
    pus.emplace_back(asym_lo.Clone());
    pus.emplace_back(asym_hi.Clone());
    pus.emplace_back(scat_anlz.Clone());
    pus.emplace_back(migrate.Clone()); // After this line, particles are all within borders.
    pus.emplace_back(export_prep.Clone());

    return pus;
  }
}

namespace particle {
  constexpr pic::real_t gamma_fd = 20.0;
  constexpr pic::real_t gamma_off = 15.0;
  constexpr pic::real_t Ndot_fd = 0.25;
  constexpr pic::real_t E_ph = 4.0;

  template < typename R >
  auto set_up() {
    map<Properties> properties;
    {
      properties.insert(species::electron, {1,-1,"electron","el"});
      //properties.insert(species::positron, {1,1,"positron","po"});
      properties.insert(species::ion, { 5, 1, "ion","io"});
      //properties.insert(species::photon, { 0, 0, "photon","ph" });
    }

    using namespace pic;
    {
      constexpr auto* lorentz = force::template lorentz<real_t,Specs,vParticle>;
      constexpr auto* landau0 = force::landau0<real_t,Specs,vParticle>;
      constexpr auto* gravity = force::gravity<real_t,Specs,vParticle>;
      real_t landau0_B_thr = 0.1;

      if ( properties.has(species::electron) ) {
        auto sp = species::electron;
        Force<real_t,Specs> force;
        const auto& prop = properties[sp];

        force.add( lorentz, ( pic::w_gyro_unitB * prop.charge_x ) / prop.mass_x );
        force.add( gravity, gravity_strength );
        // force.add( landau0, landau0_B_thr );

        force.Register(sp);
      }
      if ( properties.has(species::positron) ) {
        auto sp = species::positron;
        Force<real_t,Specs> force;
        const auto& prop = properties[sp];

        force.add( lorentz, ( pic::w_gyro_unitB * prop.charge_x ) / prop.mass_x );
        force.add( gravity, gravity_strength );
        // force.add( landau0, landau0_B_thr );

        force.Register(sp);
      }
      if ( properties.has(species::ion) ) {
        auto sp = species::ion;
        Force<real_t,Specs> force;
        const auto& prop = properties[sp];

        force.add( lorentz, ( pic::w_gyro_unitB * prop.charge_x ) / prop.mass_x );
        force.add( gravity, gravity_strength );
        // force.add( landau0, landau0_B_thr );

        force.Register(sp);
      }
    }

    using Ptc_t = typename array<real_t,Specs>::particle_type;
    {
      Scat<real_t,Specs> ep_scat;

      ep_scat.eligs.push_back([](const Ptc_t& ptc){ return ptc.q()[0] < std::log(9.0); });

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
      Scat<real_t,Specs> photon_scat;
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

#include "io/exportee.hpp"
#include "msh/mesh_shape_interplay.hpp"
#include <cassert>

namespace io {
  template < typename T, int DGrid >
  constexpr apt::array<T, DGrid> I2std ( const apt::Index<DGrid>& I ) {
    apt::array<T, DGrid> res;
    for ( int i = 0; i < DGrid; ++i )
      res[i] = I[i] + 0.5; // interpolate to MIDWAY
    return res;
  }

  template < int F, typename R, int DGrid, typename ShapeF, typename RJ >
  apt::array<R,3> field_self ( apt::Index<DGrid> I,
                               const apt::Grid<R,DGrid>& grid,
                               const field::Field<R, 3, DGrid>& E,
                               const field::Field<R, 3, DGrid>& B,
                               const field::Field<RJ, 3, DGrid>& J ) {
    if constexpr ( F == 0 ) {
        return msh::interpolate( E, I2std<R>(I), ShapeF() );
      } else if ( F == 1 ) {
      return msh::interpolate( B, I2std<R>(I), ShapeF() );
    } else if ( F == 2 ) {
      auto x = msh::interpolate( J, I2std<RJ>(I), ShapeF() );
      for ( int i = 0; i < 3; ++i ) x[i] *= field::PREJ;
      return { x[0], x[1], x[2] };
    } else {
      static_assert(F < 3);
    }
  }

  template < typename RDS, int DGrid >
  void divide_flux_by_area ( field::Field<RDS,3,DGrid>& fds, const apt::Grid<RDS,DGrid>& grid, int num_comps, const mpi::CartComm& ) {
    using Metric = metric::LogSpherical<RDS>;
    // define a function pointer.
    RDS(*hh_func)(RDS,RDS,RDS) = nullptr;
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

  template < typename R, int DGrid, typename ShapeF, typename RJ >
  apt::array<R,3> EparaB ( apt::Index<DGrid> I,
                           const apt::Grid<R,DGrid>& grid,
                           const field::Field<R, 3, DGrid>& E,
                           const field::Field<R, 3, DGrid>& B,
                           const field::Field<RJ, 3, DGrid>& J ) {
    auto B_itpl = msh::interpolate( B, I2std<R>(I), ShapeF() );
    B_itpl /= ( apt::abs(B_itpl) + 1e-16 );
    return {apt::dot( msh::interpolate( E, I2std<R>(I), ShapeF() ), B_itpl ),
            0.0, 0.0};
  }

  template < typename R, int DGrid, typename ShapeF, typename RJ >
  apt::array<R,3> EdotJ ( apt::Index<DGrid> I,
                          const apt::Grid<R,DGrid>& grid,
                          const field::Field<R, 3, DGrid>& E,
                          const field::Field<R, 3, DGrid>& B,
                          const field::Field<RJ, 3, DGrid>& J ) {
    return {apt::dot( msh::interpolate( E, I2std<R>(I), ShapeF() ),
                      msh::interpolate( J, I2std<RJ>(I), ShapeF() ) ),
            0.0, 0.0};
  }

  // Poloidal flux function, LogSpherical
  template < typename R, int DGrid, typename ShapeF, typename RJ >
  apt::array<R,3> dFlux_pol ( apt::Index<DGrid> I,
                              const apt::Grid<R,DGrid>& grid,
                              const field::Field<R, 3, DGrid>& E,
                              const field::Field<R, 3, DGrid>& B,
                              const field::Field<RJ, 3, DGrid>& J ) {
    // F = \int_Br_d\theta, we want F to be all MIDWAY. Since Br is (MIDWAY, INSITU), it is automatically the natural choice
    const auto& Br = B[0];
    return { Br(I) * std::exp( 2.0 * grid[0].absc( I[0], 0.5 ) ) * std::sin( grid[1].absc( I[1], 0.0 ) ), 0.0, 0.0 };
  }

  template < typename RDS, int DGrid >
  void integrate_dFlux ( field::Field<RDS,3,DGrid>& fds, const apt::Grid<RDS,DGrid>& grid, int num_comps, const mpi::CartComm& cart ) {
    // Flux_t - Flux_{t-1} = dFlux_t, can be in-placed
    // integrate in theta direction
    assert(num_comps == 1);
    auto dFlux = fds[0]; // TODOL semantics
    const auto& mesh = fds.mesh();
    apt::Index<DGrid> ext = apt::range::size(mesh.range());
    std::vector<RDS> buf;
    { // one value from each theta row
      std::size_t size = 1;
      for ( int i = 0; i < DGrid; ++i ) size *= ext[i];
      size /= ext[1];
      buf.reserve(size);
    }
    for ( const auto& trI : apt::project_out( 1, {}, ext ) ) {
      apt::Longidx n (1,-1);
      dFlux(trI + n) = 0.0;
      for ( ++n; n < ext[1]; ++n ) {
        dFlux(trI + n) += dFlux(trI + (n-1));
      }
      n = ext[1] - 1;
      buf.push_back(dFlux(trI + n));
    }
    // do an exclusive scan then add the scanned value back
    cart.exscan_inplace(buf.data(), buf.size());
    int idx = 0;
    for ( const auto& trI : apt::project_out( 1, {}, ext ) ) {
      auto val = buf[idx++];
      for ( apt::Longidx n (1,0); n < ext[1]; ++n ) dFlux(trI + n) += val;
    }
  }

  template < typename R, int DGrid, typename ShapeF, typename RJ >
  apt::array<R,3> pair_creation_rate ( apt::Index<DGrid> I,
                                       const apt::Grid<R,DGrid>& grid,
                                       const field::Field<R, 3, DGrid>& ,
                                       const field::Field<R, 3, DGrid>& ,
                                       const field::Field<RJ, 3, DGrid>&  ) {
    auto x = msh::interpolate( RTD<R,DGrid>::pc_counter, I2std<R>(I), ShapeF() );
    return { x[0] / RTD<R,DGrid>::pc_cumulative_time, 0, 0};
  }

  template < typename R, int DGrid, typename ShapeF, typename RJ >
  apt::array<R,3> volume_scale ( apt::Index<DGrid> I,
                                 const apt::Grid<R,DGrid>& grid,
                                 const field::Field<R, 3, DGrid>& ,
                                 const field::Field<R, 3, DGrid>& ,
                                 const field::Field<RJ, 3, DGrid>&  ) {
    return { pic::Metric::hhh(grid[0].absc( I[0],0.5), grid[1].absc( I[1],0.5) ), 0, 0 };
  }

  template < particle::species SP, typename R, int DGrid, typename ShapeF, typename RJ >
  apt::array<R,3> frac_J_sp ( apt::Index<DGrid> I,
                              const apt::Grid<R,DGrid>& grid,
                              const field::Field<R, 3, DGrid>& ,
                              const field::Field<R, 3, DGrid>& ,
                              const field::Field<RJ, 3, DGrid>& J ) {
    auto q = I2std<RJ>(I);
    auto j_sp = msh::interpolate( RTD<R,DGrid>::Jsp[SP], q, ShapeF() );
    auto j = msh::interpolate( J, q, ShapeF() );
    for ( int i = 0; i < 3; ++i ) {
      if ( j[i] == 0 ) j_sp[i] = 0;
      else j_sp[i] /= j[i];
    }
    return { j_sp[0], j_sp[1], j_sp[2] };
  }
  // void delE_v_J();

  // void delB();

  template < typename RDS,
             int DGrid,
             typename R,
             typename RJ
             >
  auto set_up_field_export() {
    std::vector<FieldExportee<RDS, DGrid, R, RJ>*> fexps;
    {
      using FA = FexpTbyFunction<RDS,DGrid,R,RJ>;
      using pic::ShapeF;

      fexps.push_back( new FA ( "E", 3,
                                field_self<0,R, DGrid, ShapeF, RJ>,
                                nullptr
                                ) );
      fexps.push_back( new FA ( "B", 3,
                                field_self<1,R, DGrid, ShapeF, RJ>,
                                nullptr
                                ) );
      fexps.push_back( new FA ( "J4X", 3,
                                field_self<2,R, DGrid, ShapeF, RJ>,
                                divide_flux_by_area<RDS, DGrid>
                                ) );
      fexps.push_back( new FA ( "EparaB", 1,
                                EparaB<R, DGrid, ShapeF, RJ>,
                                nullptr
                                ) );
      fexps.push_back( new FA ( "EdotJ", 1,
                                EdotJ<R, DGrid, ShapeF, RJ>,
                                nullptr
                                ) );
      fexps.push_back( new FA ( "Flux", 1,
                                dFlux_pol<R, DGrid, ShapeF, RJ>,
                                integrate_dFlux<RDS, DGrid>
                                ) );
      fexps.push_back( new FA ( "PairCreationRate", 1,
                                pair_creation_rate<R, DGrid, ShapeF, RJ>,
                                nullptr
                                ) );
      fexps.push_back( new FA ( "VolumeScale", 1,
                                volume_scale<R, DGrid, ShapeF, RJ>,
                                nullptr
                                ) );

      if ( RTD<R,DGrid>::is_export_Jsp ) {
        using namespace particle;
        for ( auto sp : RTD<R,DGrid>::Jsp ) {
          switch(sp) {
          case species::electron :
            fexps.push_back( new FA ( "fJ_Electron", 3,
                                      frac_J_sp<species::electron,R, DGrid, ShapeF, RJ>,
                                      divide_flux_by_area<RDS, DGrid>
                                      ) ); break;
          case species::positron :
            fexps.push_back( new FA ( "fJ_Positron", 3,
                                      frac_J_sp<species::positron,R, DGrid, ShapeF, RJ>,
                                      divide_flux_by_area<RDS, DGrid>
                                      ) ); break;
          case species::ion :
            fexps.push_back( new FA ( "fJ_Ion", 3,
                                      frac_J_sp<species::ion,R, DGrid, ShapeF, RJ>,
                                      divide_flux_by_area<RDS, DGrid>
                                      ) ); break;
          default: ;
          }
        }
      }
    }
    return fexps;
  }

}

namespace io {
  template < typename R, template < typename > class S >
  apt::array<R,3> ptc_num ( const particle::Properties& prop, const typename particle::array<R,S>::const_particle_type& ptc ) {
    return { 1.0, 0.0, 0.0 };
  }

  template < typename R, template < typename > class S >
  apt::array<R,3> ptc_energy ( const particle::Properties& prop, const typename particle::array<R,S>::const_particle_type& ptc ) {
    return { std::sqrt( (prop.mass_x != 0) + apt::sqabs(ptc.p()) ), 0.0, 0.0 };
  }

  template < typename R, template < typename > class S >
  apt::array<R,3> ptc_momentum ( const particle::Properties& prop, const typename particle::array<R,S>::const_particle_type& ptc ) {
    return { ptc.p()[0], ptc.p()[1], ptc.p()[2] };
  }

  template < int DGrid, typename RDS, template < typename > class S >
  void fold_back_at_axis ( field::Field<RDS,3,DGrid>& field, const apt::Grid<RDS,DGrid>& grid, int num_comps ) {
    // NOTE field is assumed to have all-MIDWAY offset
    for ( int i = 0; i < num_comps; ++i ) {
      for ( int dim = 0; dim < DGrid; ++dim )
        assert( field[i].offset()[dim] == MIDWAY );
    }

    constexpr int axis_dir = 1;
    constexpr auto PI = std::acos(-1.0l);
    bool is_lower = std::abs( grid[axis_dir].lower() - 0.0 ) < grid[axis_dir].delta();
    bool is_upper = std::abs( grid[axis_dir].upper() - PI ) < grid[axis_dir].delta();

    auto add_assign = []( RDS& a, RDS& b ) noexcept -> void {a += b; b = a;};
    auto sub_assign = []( RDS& a, RDS& b ) noexcept -> void {a -= b; b = -a;};

    const auto& range = field.mesh().range();
    if ( is_lower ) {
      auto rb = apt::range::begin(range);
      auto re = apt::range::end(range);
      rb[axis_dir] = -1;
      re[axis_dir] = 0; // 0 as the end is appropriate because field is MIDWAY
      axissymmetrize( field[0], add_assign, rb, re, false );
      for ( int i = 1; i < num_comps; ++i )
        axissymmetrize( field[0], sub_assign, rb, re, false );
    }
    if ( is_upper ) {
      auto rb = apt::range::begin(range);
      auto re = apt::range::end(range);
      rb[axis_dir] = re[axis_dir]; // FIXME double check. There are domain indices not global indices
      re[axis_dir] += 1;
      axissymmetrize( field[0], add_assign, rb, re, true );
      for ( int i = 1; i < num_comps; ++i )
        axissymmetrize( field[0], sub_assign, rb, re, true );
    }
  }

  template < typename RDS,
             int DGrid,
             typename R,
             template < typename > class S>
  auto set_up_particle_export() {
    constexpr int DSratio = 1;
    std::vector<PtcExportee<RDS, DGrid, R, S>*> pexps;
    {
      using PA = PexpTbyFunction<RDS,DGrid,R,S,particle::induced_shapef_t<pic::ShapeF,DSratio>>;
      pexps.push_back( new PA ("Num", 1,
                               ptc_num<R,S>,
                               fold_back_at_axis< DGrid, RDS, S >
                               ) );

      pexps.push_back( new PA ("E", 1,
                               ptc_energy<R,S>,
                               fold_back_at_axis< DGrid, RDS, S >
                               ) );

      pexps.push_back( new PA ("P", 3,
                               ptc_momentum<R,S>,
                               fold_back_at_axis< DGrid, RDS, S >
                               ) );
    }

    return pexps;
  }
}

namespace io {
  constexpr bool is_collinear_mesh = false; // FIXME this is an ad hoc fix

  template < typename R, int DGrid >
  void export_prior_hook( const apt::Grid< R, DGrid >& grid, const dye::Ensemble<DGrid>& ens ) {
    auto& pc = RTD<R,DGrid>::pc_counter;
    ens.reduce_to_chief( mpi::by::SUM, pc[0].data().data(), pc[0].data().size() );
  }

  template < typename R, int DGrid >
  void export_post_hook() {
    auto& pc = RTD<R,DGrid>::pc_counter;
    std::fill( pc[0].data().begin(), pc[0].data().end(), 0 );
    RTD<R,DGrid>::pc_cumulative_time = 0;
    if ( RTD<R,DGrid>::is_export_Jsp ) {
      // clear Jsp to save some space
      for ( auto sp : RTD<R,DGrid>::Jsp )
        RTD<R,DGrid>::Jsp[sp] = {};
    }
    RTD<R,DGrid>::skin_depth = {};
  }
}

#include <sstream>
#include "apt/print.hpp"

namespace pic {

  std::string characteristics(std::string indent) {
    std::ostringstream o;
    auto gamma_0 = std::pow(field::Omega,2.0) * w_gyro_unitB;
    o << indent << "gamma_0=" << apt::fmt("%.0f", gamma_0 ) << std::endl;
    o << indent << "w_pic dt=" << apt::fmt("%.4f", wdt_pic ) << std::endl;
    o << indent << "re=" << apt::fmt("%.4f", classic_electron_radius() ) << std::endl;
    o << indent << "Ndot_GJ=" << apt::fmt("%.4e", gamma_0 / classic_electron_radius() ) << std::endl;

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
      o << ", L_CR/L_sd=" << E_ph * Ndot_fd * field::Omega * std::pow( gamma_0 / gamma_fd, 3.0 );
    }
    return o.str();
  }
}

#endif
