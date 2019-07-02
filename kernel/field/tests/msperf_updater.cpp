#include "testfw/testfw.hpp"
#include "field/updater_impl.hpp"
#include "field/sync_impl.hpp"
#include "timer/timer.hpp"
#include "mpipp/mpi++.hpp"

TEMPLATE_TEST_CASE("Time field updater 2D", "[mpi]"
                   // periodic
                   , (aio::IndexType<-1,-1>)
                   , (aio::IndexType<-2,-2>)
                   , (aio::IndexType<-4,-4>)
                   ) {
  using Real = double;
  constexpr int DGrid = 2;
  using RealJ = double;

  constexpr auto Cartesian_Partition = TestType::get();
  CAPTURE(Cartesian_Partition);

  constexpr apt::Index<DGrid> bulk_dims = {128,128};
  constexpr mani::Grid<Real,DGrid> supergrid
    = {{ { 0.0, 1.0, bulk_dims[0] }, { 0.0, 1.0, bulk_dims[1] } }};

  auto cart_opt = aio::make_cart( Cartesian_Partition, mpi::world );
  if ( cart_opt ) {
    const auto coords = cart_opt -> coords();
    apt::array< apt::pair<bool>, DGrid > is_at_boundary;
    {
      for ( int i = 0; i < DGrid; ++i ) {
        is_at_boundary[i][LFT] = ( Cartesian_Partition[i] > 0 && 0 == coords[i] );
        is_at_boundary[i][RGT] = ( Cartesian_Partition[i] > 0 && Cartesian_Partition[i] - 1 == coords[i] );
      }
    }

    mani::Grid<Real,DGrid> grid; // local grid
    for ( int i = 0; i < DGrid; ++i )
      grid[i] = supergrid[i].divide( std::abs(Cartesian_Partition[i]), coords[i] );

    constexpr int guard = 1;
    field::Updater<Real,DGrid,RealJ> fu(*cart_opt, grid, is_at_boundary, guard);

    field::Field<Real, 3, DGrid> E {{ bulk_dims, guard }};
    field::Field<Real, 3, DGrid> B {{ bulk_dims, guard }};
    field::Field<RealJ, 3, DGrid> J {{ bulk_dims, guard }};

    // set up some random J
    aio::unif_real<RealJ> unif;
    for ( int i = 0; i < 3; ++i ) {
      for ( const auto& I : apt::Block(bulk_dims) ) {
        J[i](I) = 1000.0 * ( 2.0 * unif() - 1.0 );
      }
    }

    tmr::Timestamp t1;
    const int N = 1000;
    const Real dt = 0.005;
    for ( int i = 0; i < N; ++i ) {
      fu(E,B,J,dt,i);
    }
    auto dur = t1.lapse();
    cart_opt->barrier();
    WARN("Timing field updater per timestep = " << dur.in_units_of("ms").val() << "/" << N << " = " << dur.in_units_of("ms").val() / N << dur.unit() );
  }
  mpi::world.barrier();
}
