#ifndef _PIC_IMPL_HPP_
#define _PIC_IMPL_HPP_

#include "apt/numeric.hpp"
#include "apt/index.hpp"

#include "pic/module_range.hpp"
#include "pic/diffs/diff_cartesian.hpp"
#include "manifold/curvilinear_impl.hpp"

#include "manifold/grid.hpp"

#include "field/field.hpp"
#include "field/params.hpp"
#include "field/haugbolle_solver/updater.hpp"

#include "particle/properties.hpp"
#include "particle/map.hpp"
#include "particle/array.hpp"
#include "particle/particle.hpp"

#include "dye/ensemble.hpp"

#include "pic.hpp"

namespace pic {
  using Metric = mani::CartesianCoordSys;

  inline constexpr const char* project_name = "Cartesian256";
  inline constexpr const char* datadir_prefix = "/home/hooray/Projects/Pigeon/Data/";

  inline constexpr apt::array<int,DGrid> dims = { 1, 1 };
  inline constexpr apt::array<bool,DGrid> periodic = {true,true};
  inline constexpr real_t dt = 0.003;

  constexpr mani::Grid<real_t,DGrid> supergrid
  = {{ { 0.0, 1.0, 256 }, { 0.0, 1.0, 256 } }};

  inline constexpr real_t wdt_pic = 1.0 / 30.0;
  inline constexpr real_t w_gyro_unitB = 3750; // set the impact of unit field strength on particle

  inline void set_resume_dir( std::optional<std::string>& dir ) {}
  inline constexpr int total_timesteps = 1000;

  constexpr real_t classic_electron_radius () noexcept {
    real_t res = wdt_pic * wdt_pic / ( 4.0 * std::acos(-1.0l) * dt * dt);
    apt::foreach<0,DGrid>( [&res](const auto& g) { res *= g.delta(); }, supergrid );
    return res;
  }
}

namespace pic {
  inline constexpr ModuleRange sort_particles_mr { true, 0, 100 };

  inline constexpr ModuleRange export_data_mr { false, 0, 100 };
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
}

namespace field {
  using pic::real_t;
  constexpr int order_precision = 2; // order precision of update scheme, this means error is O(x^(order + 1));
  constexpr int order_derivative = 2; // order of derivative of update before field communication is necessary
  static_assert( order_precision % 2 == 0 );
  constexpr auto PREJ = 4.0 * std::acos(-1.0l) * pic::classic_electron_radius() / pic::w_gyro_unitB;

  constexpr int myguard = std::max(order_derivative * order_precision / 2, ( pic::ShapeF::support() + 3 ) / 2 ); // NOTE minimum number of guards of J on one side is ( supp + 3 ) / 2

  template < int DGrid, typename R, typename RJ >
  auto set_up_field_actions() {
    std::vector<std::unique_ptr<Action<R,DGrid,RJ>>> fus;

    Haugbolle<R,DGrid,RJ> fu_bulk;
    {
      auto& fu = fu_bulk;
      fu.setIb( { 0, 0 } );
      fu.setIe({ pic::supergrid[0].dim(), pic::supergrid[1].dim() });
      fu.setGuard ({ myguard, myguard, myguard, myguard });
      fu.set_preJ(PREJ);

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
    struct InitialCondition : public apt::ActionBase<DGrid> {
      InitialCondition* Clone() const override { return new InitialCondition(*this); }

      void operator() ( const mani::Grid<R,DGrid>& grid,
                        field::Field<R, 3, DGrid>& E,
                        field::Field<R, 3, DGrid>& B,
                        field::Field<RJ, 3, DGrid>& J,
                        particle::map<particle::array<R,S>>& particles
                        ) const {
      }

      int initial_timestep() const noexcept { return 0; }
    } ic;
    ic.setIb({0,0});
    ic.setIe({supergrid[0].dim(), supergrid[1].dim()+1});

    return ic;
  }
}

namespace particle {
  template < int DGrid, typename R, template < typename > class S,
             typename ShapeF, typename RJ >
  auto set_up_particle_actions() {
    std::vector<std::unique_ptr<Action<DGrid,R,S,RJ>>> pus;

    Updater<DGrid,R,S,ShapeF,RJ> pu;
    pu.set_update_q(pic::Metric::geodesic_move<apt::vVec<R,3>, apt::vVec<R,3>, R>);

    pus.emplace_back(pu.Clone());
    return pus;
  }
}

namespace particle {
  template < typename R = double > // TODOL this is an ad hoc trick to prevent calling properties in routines where it is not needed
  void set_up_properties() {
    properties.insert(species::electron, {1,-1,"electron","el"});
    properties.insert(species::positron, {1,1,"positron","po"});
  }

  // NOTE called in particle updater
  template < typename R >
  void set_up() {
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
      if ( properties.has(species::positron) ) {
        auto sp = species::positron;
        Force<real_t,Specs> force;
        const auto& prop = properties[sp];

        force.add( lorentz, ( pic::w_gyro_unitB * prop.charge_x ) / prop.mass_x );

        force.Register(sp);
      }
    }

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
                               const mani::Grid<R,DGrid>& grid,
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

  template < typename RDS, typename ShapeF, int DGrid, typename R, typename RJ>
  auto set_up_field_export() {
    std::vector<FieldExportee<RDS, DGrid, R, ShapeF, RJ>*> fexps;
    {
      using FA = FieldAction<RDS,DGrid,R,ShapeF,RJ>;

      fexps.push_back( new FA ( "E", 3, field_self<0,R, DGrid, ShapeF, RJ>, nullptr) );
      fexps.push_back( new FA ( "B", 3, field_self<1,R, DGrid, ShapeF, RJ>, nullptr) );
      fexps.push_back( new FA ( "J4X", 3, field_self<2,R, DGrid, ShapeF, RJ>, nullptr) );
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

  template < typename RDS,
             typename ShapeF,
             int DGrid,
             typename R,
             template < typename > class S>
  auto set_up_particle_export() {
    std::vector<PtcExportee<RDS, DGrid, R, S, ShapeF>*> pexps;
    {
      using PA = PtcAction<RDS,DGrid,R,S,ShapeF>;
      pexps.push_back( new PA ("Num", 1, ptc_num<R,S>, nullptr) );
      pexps.push_back( new PA ("E", 1, ptc_energy<R,S>, nullptr ) );
      pexps.push_back( new PA ("P", 3, ptc_momentum<R,S>, nullptr ) );
    }

    return pexps;
  }
}

namespace io {
  constexpr bool is_collinear_mesh = true; // TODO this is an ad hoc fix
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
