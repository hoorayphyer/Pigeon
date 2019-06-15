#include "testfw/testfw.hpp"
#include "particle/updater_impl.hpp"
#include "timer/timer.hpp"
#include "mpipp/mpi++.hpp"
#include "gen.hpp"

using namespace particle;
using namespace pic;
SCENARIO("Time particle updater", "[particle]") {
  util::Rng<real_t> rng; // TODO seed
  using PU_t = Updater<DGrid, real_t, Specs, ShapeF, real_j_t, Metric>;

  constexpr apt::Index<DGrid> bulk_dims = {128,128};
  constexpr mani::Grid<Real,DGrid> grid
    = {{ { 0.0, 1.0, bulk_dims[0] }, { 0.0, 1.0, bulk_dims[1] } }};

  PU_t pu( grid, rng );

  tmr::Timestamp t1;
  const int N = 1000;
  const Real dt = 0.005;
  for ( int i = 0; i < N; ++i ) {
    // fu(E,B,J,dt,i);
  }
  auto dur = t1.lapse();
  WARN("Timing particle updater per timestep = " << dur << "/" << N << " = " << dur / N << " ms" );
}
