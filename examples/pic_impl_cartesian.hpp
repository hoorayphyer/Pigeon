#ifndef _PIC_IMPL_HPP_
#define _PIC_IMPL_HPP_

#include "../examples/pic_prior_impl.hpp"
#include "apt/numeric.hpp"
#include "pic/module_range.hpp"
#include "metric/cartesian.hpp"
#include "io/exportee_by_function.hpp"


namespace pic {
  using Metric = metric::Cartesian<real_t>;

  inline constexpr const char* project_name = "Cartesian";
  inline constexpr const char* datadir_prefix = "../Data/";

  inline constexpr apt::array<int,DGrid> dims = { 1, 1 };
  inline constexpr apt::array<bool,DGrid> periodic = {true,true};
  inline constexpr real_t dt = 0.01;

  constexpr Grid supergrid
  = {{ { 0.0, 1.0, 64 }, { 0.0, 1.0, 64 } }};

  inline constexpr real_t wdt_pic = 1.0 / 30.0;
  inline constexpr real_t w_gyro_unitB = 3750; // set the impact of unit field strength on particle

  inline void set_resume_dir( std::optional<std::string>& dir ) {}
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

namespace pic {
  constexpr int order_precision = 2; // order precision of update scheme, this means error is O(x^(order + 1));
  constexpr int number_iteration = 4; // number of iterations in inverting the operator in Haugbolle
  static_assert( order_precision % 2 == 0 );
  constexpr auto PREJ = 4.0 * std::acos(-1.0l) * pic::classic_electron_radius() / pic::w_gyro_unitB;

  constexpr int myguard = std::max(1 + order_precision * (1+number_iteration) / 2, ( pic::ShapeF::support() + 3 ) / 2 ); // NOTE minimum number of guards of J on one side is ( supp + 3 ) / 2
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

    map<unsigned int> trace_counter {};

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
  namespace yee = ::field::yee;
  template < ::field::offset_t Ftype >
  constexpr auto Diff = ::field::Diff<DGrid,real_t,Ftype>;
  constexpr auto diff_one = ::field::diff_one<real_t>;

  auto set_up_field_actions() {
    std::vector<std::unique_ptr<FieldAction>> fus;
    namespace range = apt::range;

    Haugbolle fu_bulk;
    {
      auto& fu = fu_bulk;
      fu.setName("Bulk");
      fu[0] = { 0, pic::supergrid[0].dim(), myguard };
      fu[1] = { 0, pic::supergrid[1].dim(), myguard };
      fu.set_fourpi(PREJ);
      fu.set_number_iteration(number_iteration);
      fu.set_implicit(0.8); // 0.8 is stable for E_phi FIXME not really, grows after t = 15

      for ( int i = 0; i < 3; ++i ) { // i is coordinate
        for ( int j = 0; j < 3; ++j ) { // j is Fcomp
          if ( i == j ) continue;
          fu.set_D(yee::Etype,j,i, Diff<yee::Etype>(j,i) );
          fu.set_D(yee::Btype,j,i, Diff<yee::Btype>(j,i) );
        }
      }
      for ( int Ftype = Ftype; Ftype < 2; ++Ftype )
        for ( int i = 0; i < 3; ++i )
          fu.set_hh( Ftype,i, diff_one );
    }

    fus.emplace_back(fu_bulk.Clone());

    return fus;
  }
}

namespace pic {
  auto set_up_particle_properties() {
    map<Properties> properties;
    {
      properties.insert(species::electron, {1,-1,"electron","el"});
      properties.insert(species::ion, { 5, 1, "ion","io"});
    }

    {
      constexpr auto* lorentz = ::particle::force::template lorentz<real_t,Specs,::particle::vParticle>;

      if ( properties.has(species::electron) ) {
        auto sp = species::electron;
        Force force;
        const auto& prop = properties[sp];

        force.add( lorentz, ( w_gyro_unitB * prop.charge_x ) / prop.mass_x );

        force.Register(sp);
      }
      if ( properties.has(species::ion) ) {
        auto sp = species::ion;
        Force force;
        const auto& prop = properties[sp];

        force.add( lorentz, ( pic::w_gyro_unitB * prop.charge_x ) / prop.mass_x );

        force.Register(sp);
      }
    }

    return properties;
  }

}

namespace pic {
  using R = real_t;

  auto set_up_particle_actions() {
    namespace range = apt::range;
    std::vector<std::unique_ptr<PtcAction>> pus;

    PtcUpdater pu;
    {
      pu.setName("MainUpdate");
      pu.set_update_q(Metric::geodesic_move<apt::vVec<R,3>, apt::vVec<R,3>>);
    }

    ::particle::Migrator<DGrid,real_t,Specs,ShapeF,real_j_t> migrate;
    {
      migrate.setName("MigrateParticles");
      migrate.set_supergrid(pic::supergrid);
    }

    pus.emplace_back(pu.Clone());
    pus.emplace_back(migrate.Clone());

    return pus;
  }
}

namespace pic {
  auto set_up_initial_conditions() {
    // local class in a function
    struct InitialCondition : public apt::ActionBase<DGrid> {
      InitialCondition* Clone() const override { return new InitialCondition(*this); }

      void operator() ( const Grid& grid, Field<3>& , Field<3>& B,
                        JField& , map<PtcArray>& particles) const {
        // particles[species::electron].emplace_back({0.5,0.5,0.0}, {10.0,0.0,0.0}, 1.0, species::electron);
      }

    } ic;
    ic[0] = { 0, supergrid[0].dim() };
    ic[1] = { 0, supergrid[1].dim() };

    return ic;
  }
}

#include "io/exportee.hpp"
#include "msh/mesh_shape_interplay.hpp"
#include <cassert>

namespace pic {
  constexpr bool is_collinear_mesh = true; // FIXME this is an ad hoc fix

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
      auto x = msh::interpolate( J, I2std(I), ShapeF() );
      for ( int i = 0; i < 3; ++i ) x[i] *= PREJ;
      return { x[0], x[1], x[2] };
    } else {
      static_assert(F < 3);
    }
  }

  apt::array<real_t,3> EparaB ( Index I, const Grid& grid, const Field<3>& E,
                                const Field<3>& B, const JField& J ) {
    auto B_itpl = msh::interpolate( B, I2std(I), ShapeF() );
    B_itpl /= ( apt::abs(B_itpl) + 1e-16 );
    return {apt::dot( msh::interpolate( E, I2std(I), ShapeF() ), B_itpl ),
            0.0, 0.0};
  }

  apt::array<real_t,3> EdotJ ( Index I, const Grid& grid, const Field<3>& E,
                               const Field<3>& B, const JField& J ) {
    return {apt::dot( msh::interpolate( E, I2std(I), ShapeF() ),
                      msh::interpolate( J, I2std(I), ShapeF() ) ),
            0.0, 0.0};
  }


  auto set_up_field_export() {
    std::vector<::io::FieldExportee<real_export_t, DGrid, real_t, real_j_t>*> fexps;
    {
      using FA = ::io::FexpTbyFunction<real_export_t, DGrid, real_t, real_j_t>;

      fexps.push_back( new FA ( "E", 3, field_self<0>, nullptr) );
      fexps.push_back( new FA ( "B", 3, field_self<1>, nullptr) );
      fexps.push_back( new FA ( "J4X", 3, field_self<2>, nullptr) );
      fexps.push_back( new FA ( "EparaB", 1, EparaB, nullptr) );
      fexps.push_back( new FA ( "EdotJ", 1, EdotJ, nullptr) );
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
  void export_post_hook() {}
}

#include <sstream>
#include "apt/print.hpp"

namespace pic {

  std::string characteristics(std::string indent) {
    std::ostringstream o;
    return o.str();
  }
}

#endif
