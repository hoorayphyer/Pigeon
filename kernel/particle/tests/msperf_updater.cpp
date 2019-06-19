#include "testfw/testfw.hpp"
#include "manifold/curvilinear.hpp"
#include "particle/shapef.hpp"
#include "particle/array_impl.hpp"
#include "particle/forces_impl.hpp"
#include "particle/scattering_impl.hpp"
#include "msh/mesh_shape_interplay_impl.hpp"
#include "msh/current_deposition_impl.hpp"
#include "particle/updater_impl.hpp"
#include "timer/timer.hpp"
#include "particle/properties.hpp"

using namespace particle;

namespace particle {
  map<Properties> properties;
}

// TODOL users may forget to sync value and state. Add another layer then
template < typename T >
struct PtcSpecs {
  using value_type = T;
  static constexpr int Dim = 3;
  using state_type = apt::copy_cvref_t<T,unsigned long long>;
};

SCENARIO("Time particle updater", "[particle]") {
  constexpr int DGrid = 2;
  using Real = double;
  using Real_J = double;
  using Metric = mani::coord<mani::coordsys::LogSpherical>;
  using ShapeF = particle::shapef_t<particle::shape::Cloud_In_Cell>;

  properties[species::electron] = {1,-1,"electron"};
  // properties[species::positron] = {1,1,"positron"};

  using Force = force::Force<Real,PtcSpecs,vParticle>;
  ForceGen<Real,PtcSpecs,vParticle> fgen;
  constexpr auto* lorentz = force::template lorentz<Real,PtcSpecs,vParticle>;
  // constexpr auto* landau0 = force::landau0<real_t,Specs,vParticle>;

  // constexpr Real gravity_strength = 1.8;
  // constexpr auto* gravity = force::gravity<real_t,Specs,vParticle>;
  {
    auto sp = species::electron;
    Force force;
    const auto& prop = properties.at(sp);

    force.add( lorentz, static_cast<Real>(prop.charge_x) / prop.mass_x );
    // force.add( gravity, gravity_strength );
    // force.add( landau0, landau0_B_thr );

    fgen.Register( sp, force );
  }
  // {
  //   auto sp = species::positron;
  //   Force force;
  //   const auto& prop = properties.at(sp);

  //   force.add( lorentz, static_cast<Real>(prop.charge_x) / prop.mass_x );
  //   // force.add( gravity, gravity_strength );
  //   // force.add( landau0, landau0_B_thr );

  //   fgen.Register( sp, force );
  // }

  util::Rng<Real> rng; // TODO seed
  using PU_t = Updater<DGrid, Real, PtcSpecs, ShapeF, Real_J, Metric>;

  constexpr apt::Index<DGrid> bulk_dims = {128,128};
  constexpr mani::Grid<Real,DGrid> grid
    = {{ { 0.0, 1.0, bulk_dims[0] }, { 0.0, 1.0, bulk_dims[1] } }};
  constexpr int guard = 1;

  PU_t pu( grid, rng, properties, fgen, {} );

  field::Field<Real, 3, DGrid> E;
  field::Field<Real, 3, DGrid> B;
  field::Field<Real_J, 3, DGrid> J;
  particle::map<particle::array<Real, PtcSpecs>> particles;

  {
    auto all_but =
      []( field::offset_t all, int but_comp ) noexcept
      {
        apt::array<field::offset_t, DGrid> res;
        apt::foreach<0,DGrid>( [&all]( auto& ofs ) { ofs = all; }, res );
        if ( but_comp < DGrid ) res[but_comp] = !all;
        return res;
      };
    E = {{ bulk_dims, guard }};
    B = {{ bulk_dims, guard }};
    // NOTE minimum number of guards of J on one side is ( supp + 1 ) / 2 + 1
    J = {{ bulk_dims, ( ShapeF::support() + 3 ) / 2 }};

    for( int i = 0; i < 3; ++i ) {
      E.set_offset( i, all_but( MIDWAY, i ) );
      B.set_offset( i, all_but( INSITU, i ) );
      J.set_offset( i, all_but( MIDWAY, i ) );
    }
    E.reset();
    B.reset();
    J.reset();
  }


  { // measure
    const Real dt = 0.005;

    tmr::Timestamp t1;
    const int Nptc = 10000;
    {
      auto& ptcs = particles[species::electron];
      for ( int i = 0; i < Nptc; ++i ) {
        ptcs.push_back(Particle<Real,PtcSpecs>(apt::Vec<Real,PtcSpecs<Real>::Dim>{0.5, 0.5, 0.5}, apt::Vec<Real,PtcSpecs<Real>::Dim>{10.0, 10.0, 10.0}, species::electron));
      }
    }

    const int N = 100;
    for ( int i = 0; i < N; ++i ) {
      pu( particles, J, E, B, dt, i );
    }
    auto dur = t1.lapse();
    dur.in_units_of("ms");
    WARN("Timing particle update per timestep per million particle = " << dur.val() * 1e6 / ( N * Nptc ) << " " << dur.unit() );
  }
}
