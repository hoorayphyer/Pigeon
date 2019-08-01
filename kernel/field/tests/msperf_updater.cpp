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
  using Real = float;
  constexpr int DGrid = 2;
  using RealJ = float;

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
    constexpr auto omega_t = [] ( Real t ) -> Real { return 0.2 * ( std::min<Real>( t / 4.0, 1.0) ); };
    constexpr Real re_over_w_gyro_unitB = 1.0;
    {
      using namespace field::ofs;
      magnetic_pole = 2; // 1 for mono-, 2 for di-
      indent = { 5, 43, guard, guard };
      damping_rate = 10.0;
    }

    field::Updater<Real,DGrid,RealJ> fu(*cart_opt, grid, is_at_boundary, guard, omega_t, re_over_w_gyro_unitB);

    field::Field<Real, 3, DGrid> E {{ bulk_dims, guard }};
    field::Field<Real, 3, DGrid> B {{ bulk_dims, guard }};
    field::Field<RealJ, 3, DGrid> J {{ bulk_dims, guard }};

    // set up some random J
    aio::unif_real<RealJ> unif;
    for ( int i = 0; i < 3; ++i ) {
      for ( const auto& I : apt::Block(bulk_dims) ) {
        J[i](I) = ( 2.0 * unif() - 1.0 ) * 0.1;
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
