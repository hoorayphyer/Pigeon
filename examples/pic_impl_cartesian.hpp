#ifndef _PIC_IMPL_HPP_
#define _PIC_IMPL_HPP_

#include "apt/numeric.hpp"
#include "apt/index.hpp"
#include "apt/grid.hpp"

#include "pic/module_range.hpp"
#include "metric/cartesian.hpp"

#include "field/field.hpp"

#include "field/haugbolle_solver/updater.hpp"

#include "particle/properties.hpp"
#include "particle/map.hpp"
#include "particle/array.hpp"
#include "particle/particle.hpp"

#include "dye/ensemble.hpp"

#include "io/exportee_by_function.hpp"

#include "pic.hpp"

namespace pic {
  using Metric = metric::Cartesian<real_t>;

  inline constexpr const char* project_name = "Cartesian";
  inline constexpr const char* datadir_prefix = "../Data/";

  inline constexpr apt::array<int,DGrid> dims = { 1, 1 };
  inline constexpr apt::array<bool,DGrid> periodic = {true,true};
  inline constexpr real_t dt = 0.01;

  constexpr apt::Grid<real_t,DGrid> supergrid
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

namespace field {
  using pic::real_t;

  constexpr int order_precision = 2; // order precision of update scheme, this means error is O(x^(order + 1));
  constexpr int number_iteration = 4; // number of iterations in inverting the operator in Haugbolle
  static_assert( order_precision % 2 == 0 );
  constexpr auto PREJ = 4.0 * std::acos(-1.0l) * pic::classic_electron_radius() / pic::w_gyro_unitB;

  // FIXME NOTE the 1 + 
  constexpr int myguard = std::max(1 + order_precision * (1+number_iteration) / 2, ( pic::ShapeF::support() + 3 ) / 2 ); // NOTE minimum number of guards of J on one side is ( supp + 3 ) / 2
}

namespace field {
  template < int DGrid, typename R, typename RJ >
  auto set_up_field_actions() {
    std::vector<std::unique_ptr<Action<R,DGrid,RJ>>> fus;
    namespace range = apt::range;

    Haugbolle<R,DGrid,RJ> fu_bulk;
    {
      auto& fu = fu_bulk;
      fu.setName("Bulk");
      fu[0] = { 0, pic::supergrid[0].dim(), myguard };
      fu[1] = { 0, pic::supergrid[1].dim(), myguard };
      fu.set_fourpi(PREJ);
      fu.set_number_iteration(number_iteration);
      fu.set_implicit(0.8);

      for ( int i = 0; i < 3; ++i ) { // i is coordinate
        for ( int j = 0; j < 3; ++j ) { // j is Fcomp
          if ( i == j ) continue;
          fu.set_D(yee::Etype,j,i, Diff<DGrid,R,yee::Etype>(j,i) );
          fu.set_D(yee::Btype,j,i, Diff<DGrid,R,yee::Btype>(j,i) );
        }
      }
      for ( int Ftype = Ftype; Ftype < 2; ++Ftype )
        for ( int i = 0; i < 3; ++i )
          fu.set_hh( Ftype,i, diff_one<DGrid,R> );
    }

    fus.emplace_back(fu_bulk.Clone());

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
                        field::Field<R, 3, DGrid>& E,
                        field::Field<R, 3, DGrid>& B,
                        field::Field<RJ, 3, DGrid>& J,
                        particle::map<particle::array<R,S>>& particles
                        ) const {
        using namespace particle;
        // particles[species::electron].emplace_back({0.5,0.5,0.0}, {10.0,0.0,0.0}, 1.0, species::electron);
      }

    } ic;
    ic[0] = { 0, supergrid[0].dim() };
    ic[1] = { 0, supergrid[1].dim() };

    return ic;
  }
}

namespace {
  template < typename R, int D >
  struct RTD { // runtime data
    static particle::map<R> N_scat;
    static field::Field<R,1,D> pc_counter;
    static R pc_cumulative_time;

    static particle::map<unsigned int> trace_counter;

    static void init( const particle::map<particle::Properties>& properties, const apt::Grid< R, D >& localgrid ) {
      for ( auto sp : properties )
        N_scat.insert( sp, 0 );

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
}

namespace particle {
  template < int DGrid, typename R, template < typename > class S,
             typename ShapeF, typename RJ >
  auto set_up_particle_actions() {
    namespace range = apt::range;
    std::vector<std::unique_ptr<Action<DGrid,R,S,RJ>>> pus;

    Updater<DGrid,R,S,ShapeF,RJ> pu;
    {
      pu.setName("MainUpdate");
      pu.set_update_q(pic::Metric::geodesic_move<apt::vVec<R,3>, apt::vVec<R,3>>);
    }

    Migrator<DGrid,R,S,ShapeF,RJ> migrate;
    {
      migrate.setName("MigrateParticles");
      migrate.set_supergrid(pic::supergrid);
    }

    pus.emplace_back(pu.Clone());
    pus.emplace_back(migrate.Clone());

    return pus;
  }
}


namespace particle {
  template < typename R >
  auto set_up() {
    map<Properties> properties;
    {
      properties.insert(species::electron, {1,-1,"electron","el"});
      properties.insert(species::ion, { 5, 1, "ion","io"});
    }

    using namespace pic;
    {
      constexpr auto* lorentz = force::template lorentz<real_t,Specs,vParticle>;
      if ( properties.has(species::electron) ) {
        auto sp = species::electron;
        Force<real_t,Specs> force;
        const auto& prop = properties[sp];

        force.add( lorentz, ( pic::w_gyro_unitB * prop.charge_x ) / prop.mass_x );

        force.Register(sp);
      }
      if ( properties.has(species::ion) ) {
        auto sp = species::ion;
        Force<real_t,Specs> force;
        const auto& prop = properties[sp];

        force.add( lorentz, ( pic::w_gyro_unitB * prop.charge_x ) / prop.mass_x );

        force.Register(sp);
      }
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
                                nullptr
                                ) );
      fexps.push_back( new FA ( "EparaB", 1,
                                EparaB<R, DGrid, ShapeF, RJ>,
                                nullptr
                                ) );
      fexps.push_back( new FA ( "EdotJ", 1,
                                EdotJ<R, DGrid, ShapeF, RJ>,
                                nullptr
                                ) );
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

  template < typename RDS, int DGrid, typename R, template < typename > class S>
  auto set_up_particle_export() {
    constexpr int DSratio = 1;
    std::vector<PtcExportee<RDS, DGrid, R, S>*> pexps;
    {
      using PA = PexpTbyFunction<RDS,DGrid,R,S,particle::induced_shapef_t<pic::ShapeF,DSratio>>;
      pexps.push_back( new PA ("Num", 1, ptc_num<R,S>, nullptr) );

      pexps.push_back( new PA ("E", 1, ptc_energy<R,S>, nullptr) );

      pexps.push_back( new PA ("P", 3, ptc_momentum<R,S>, nullptr) );
    }

    return pexps;
  }
}

namespace io {
  constexpr bool is_collinear_mesh = true; // FIXME this is an ad hoc fix

  template < typename R, int DGrid >
  void export_prior_hook( const apt::Grid< R, DGrid >& grid, const dye::Ensemble<DGrid>& ens ) {}

  template < typename R, int DGrid >
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
