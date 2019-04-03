#include "all_in_one.hpp"
#include "field/communication.cpp"
#include "parallel/mpi++.hpp"
#include <algorithm> // for std::min

// NOTE Notation: XxYxZ is the cartesian partition. Each one of X,Y,Z can be positive ( meaning nonperiodic ) or negative ( meaning periodic ).
// RATIONALE each node has a field::Field<int,3,2>, whose bulk is filled with its linearized cartesian carcoordinate.
// TODOL test on DGrid = 3
using namespace field;
using VF = Field<double,3,2>;

// NOTE this function should be nondegenerate with respect to the inputs
double field_value( const std::vector<int>& coords, int comp ) {
  double res = 0;
  for ( int i = 0; i < coords.size(); ++i )
    res += coords[i] * std::exp( ( comp + 1.0 ) / ( i + 1.0 ) );
  return res;
}

auto make_cart( std::vector<int> dims, std::vector<bool> periodic ) {
  std::optional<mpi::CartComm> cart;
  int size = 1;
  for ( auto x : dims ) size *= x;
  if ( mpi::world.size() >= size ) {
    auto comm = mpi::world.split( (mpi::world.rank() < size) );
    if ( mpi::world.rank() < size )
      cart.emplace( *comm, dims, periodic );
  }

  return cart;
}

void test_sync_guard ( const Mesh<2>& mesh, const mpi::CartComm& cart ) {
  const auto my_coords = cart.coords();
  VF vf( mesh );

  // only bulk has nonzero values
  for ( int i = 0; i < VF::NDim; ++i ) {
    double val = field_value( my_coords, i );
    for( auto I_bulk : apt::Block<2>({ mesh.bulk_dim(0), mesh.bulk_dim(1) }) ) {
      vf[i](I_bulk) = val;
    }
  }

  sync_guard_cells_from_bulk(vf, cart);

  // return Ib_bulk and ext
  auto region_begin_extent =
    []( int region_coord, int mesh_ext, int g ) -> std::tuple<int,int>
    {
     switch (region_coord) {
     case -1 : return {-g, g};
     case 0 : return {0, mesh_ext - 2 * g};
     case 1 : return {mesh_ext - 2*g, g};
     }
    };

  apt::Index<2> Ib;
  apt::Index<2> ext;
  apt::Index<2> region;
  // region[ith_dim] can be -1, 0, 1, corresponding to in left guard, in bulk, in right guard
  for ( region[1] = -1; region[1] < 2; ++region[1] ) {
    std::tie(Ib[1], ext[1]) = region_begin_extent( region[1], mesh.extent()[1], mesh.guard() );
    for ( region[0] = -1; region[0] < 2; ++region[0] ) {
      std::tie(Ib[0], ext[0]) = region_begin_extent( region[0], mesh.extent()[0], mesh.guard() );
      auto neigh = my_coords;
      neigh[0] += region[0];
      neigh[1] += region[1];

      bool at_bdry = false;
      for ( int i = 0; i < 2; ++i ) {
        if ( 0 == region[i] ) continue;
        at_bdry = ( at_bdry || !(cart.shift(i)[ 1 == region[i] ]) );
      }

      for ( int i_fld_dim = 0; i_fld_dim < VF::NDim; ++i_fld_dim ) {
        double region_val = at_bdry ? 0.0 : field_value( neigh, i_fld_dim );
        CAPTURE( mpi::world.rank(), i_fld_dim, Ib, ext, region, at_bdry );
        for ( auto I : apt::Block( ext ) )
          REQUIRE( vf[i_fld_dim](I + Ib) == region_val );
      }
    }
  }
}

void test_merge_guard ( const Mesh<2>& mesh, const mpi::CartComm& cart ) {
  const auto my_coords = cart.coords();
  VF vf( mesh );

  // all cells have nonzero values
  for ( int i = 0; i < VF::NDim; ++i ) {
    double val = field_value( my_coords, i );
    for( auto& x : vf[i].data() ) x = val;
  }

  merge_guard_cells_into_bulk(vf, cart);

  { // Check bulk value

    // return Ib_bulk and ext
    auto region_begin_extent =
      []( int region_coord, int mesh_ext, int g ) -> std::tuple<int,int>
      {
       switch (region_coord) {
       case -1 : return {0, g};
       case 0 : return {g, mesh_ext - 4 * g};
       case 1 : return {mesh_ext - 3*g, g};
       }
      };

    apt::Index<2> Ib;
    apt::Index<2> ext;
    apt::Index<2> region;
    // region[ith_dim] can be -1, 0, 1, corresponding to left part in bulk of width guard, bulk, right part in bulk of width guard
    for ( region[1] = -1; region[1] < 2; ++region[1] ) {
      std::tie(Ib[1], ext[1]) = region_begin_extent( region[1], mesh.extent()[1], mesh.guard() );
      for ( region[0] = -1; region[0] < 2; ++region[0] ) {
        std::tie(Ib[0], ext[0]) = region_begin_extent( region[0], mesh.extent()[0], mesh.guard() );


        for ( int i_fld_dim = 0; i_fld_dim < VF::NDim; ++i_fld_dim ) {
          double region_val = 0.0;
          // RATIONALE for (regx, regy), contributing neighbors are set of { ( regx, 0 ), ( regx, regy ), (0, regy), (0,0) } with duplicates all dropped. Note (0,0) is self. When reg_i = 1 or -1 corresponds to topology boundary, that dimension is ignored.
          for ( int j = std::min(region[1],0); j <= std::max(region[1],0); ++j ) {
            if ( j != 0 && !cart.shift(1)[1 == j] ) continue; // check if at cartesian topology end
            for ( int i = std::min(region[0],0); i <= std::max(region[0],0); ++i ) {
              if ( i != 0 && !cart.shift(0)[1 == i] ) continue; // check if at cartesian topology end
              auto neigh = my_coords;
              neigh[0] += i;
              neigh[1] += j;
              region_val += field_value( neigh, i_fld_dim );
            }
          }

          INFO("checking bulk areas in merge_guard");
          CAPTURE( mpi::world.rank(), i_fld_dim, Ib, ext, region  );
          // NOTE since there is addition involved, we have to use Approx
          for ( auto I : apt::Block( ext ) )
            REQUIRE( vf[i_fld_dim](I + Ib) == Approx(region_val) );
        }
      }
    }
  }

  { // Check guard cells as they need to be all cleared to avoid double depositing. // NOTE this check is very much the same as that in test_sync_guard
    // return Ib_bulk and ext
    auto region_begin_extent =
      []( int region_coord, int mesh_ext, int g ) -> std::tuple<int,int>
      {
       switch (region_coord) {
       case -1 : return {-g, g};
       case 0 : return {0, mesh_ext - 2 * g};
       case 1 : return {mesh_ext - 2*g, g};
       }
      };

    apt::Index<2> Ib;
    apt::Index<2> ext;
    apt::Index<2> region;
    // region[ith_dim] can be -1, 0, 1, corresponding to left part in bulk of width guard, bulk, right part in bulk of width guard
    for ( region[1] = -1; region[1] < 2; ++region[1] ) {
      std::tie(Ib[1], ext[1]) = region_begin_extent( region[1], mesh.extent()[1], mesh.guard() );
      for ( region[0] = -1; region[0] < 2; ++region[0] ) {
        std::tie(Ib[0], ext[0]) = region_begin_extent( region[0], mesh.extent()[0], mesh.guard() );
        if ( region[0] == 0 && region[1] == 0 ) continue; // only check guard cells
        bool at_bdry = false;
        for ( int i = 0; i < 2; ++i ) {
          if ( 0 == region[i] ) continue;
          at_bdry = ( at_bdry || !(cart.shift(i)[ 1 == region[i] ]) );
        }
        // NOTE when the guard cell region is at the cartesian topology boundary. The values are not checked. This is because these values in practice are first taken care of by boundary conditions, which will set them to zero.
        if ( at_bdry ) continue;

        for ( int i_fld_dim = 0; i_fld_dim < VF::NDim; ++i_fld_dim ) {
          INFO("checking guard cells in merge_guard");
          CAPTURE( mpi::world.rank(), i_fld_dim, Ib, ext, region, at_bdry );
          for ( auto I : apt::Block( ext ) )
            REQUIRE( vf[i_fld_dim](I + Ib) == 0.0 );
        }
      }
    }
  }
}


SCENARIO("test sync guard on nonperiodic cartesian topology", "[field][mpi]") {
constexpr auto mesh = Mesh<2>( {4, 4}, 1 );
  WHEN("1x1") {
    auto cart_opt = make_cart( {1, 1}, {false, false} );
    if ( cart_opt ) test_sync_guard(mesh, *cart_opt);
    mpi::world.barrier();
  }

  WHEN("2x1") {
    auto cart_opt = make_cart( {2, 1}, {false, false} );
    if ( cart_opt ) test_sync_guard(mesh, *cart_opt);
    mpi::world.barrier();
  }

  WHEN("1x2") {
    auto cart_opt = make_cart( {1, 2}, {false, false} );
    if ( cart_opt ) test_sync_guard(mesh, *cart_opt);
    mpi::world.barrier();
  }

  WHEN("2x2") {
    auto cart_opt = make_cart( {2, 2}, {false, false} );
    if ( cart_opt ) test_sync_guard(mesh, *cart_opt);
    mpi::world.barrier();
  }

  WHEN("4x4") {
    auto cart_opt = make_cart( {4, 4}, {false, false} );
    if ( cart_opt ) test_sync_guard(mesh, *cart_opt);
    mpi::world.barrier();
  }

  WHEN("8x8") {
    auto cart_opt = make_cart( {8, 8}, {false, false} );
    if ( cart_opt ) test_sync_guard(mesh, *cart_opt);
    mpi::world.barrier();
  }

  WHEN("16x16") {
    auto cart_opt = make_cart( {16, 16}, {false, false} );
    if ( cart_opt ) test_sync_guard(mesh, *cart_opt);
    mpi::world.barrier();
  }
}

SCENARIO("test sync guard on cartesian topology involving periodic boundary", "[field][mpi]") {
  // TODOL
}

SCENARIO("test merge guard on nonperiodic cartesian topology", "[field][mpi]") {
  constexpr auto mesh = Mesh<2>( {4, 4}, 1 );
  THEN("left and right deposited areas should not overlap") {
    for ( int i = 0; i < 2; ++i )
      REQUIRE( mesh.extent()[i] >= 4 * mesh.guard() );
  }
  WHEN("1x1") {
    auto cart_opt = make_cart( {1, 1}, {false, false} );
    if ( cart_opt ) test_merge_guard(mesh, *cart_opt);
    mpi::world.barrier();
  }

  WHEN("2x1") {
    auto cart_opt = make_cart( {2, 1}, {false, false} );
    if ( cart_opt ) test_merge_guard(mesh, *cart_opt);
    mpi::world.barrier();
  }

  WHEN("1x2") {
    auto cart_opt = make_cart( {1, 2}, {false, false} );
    if ( cart_opt ) test_merge_guard(mesh, *cart_opt);
    mpi::world.barrier();
  }

  WHEN("2x2") {
    auto cart_opt = make_cart( {2, 2}, {false, false} );
    if ( cart_opt ) test_merge_guard(mesh, *cart_opt);
    mpi::world.barrier();
  }

  WHEN("4x4") {
    auto cart_opt = make_cart( {4, 4}, {false, false} );
    if ( cart_opt ) test_merge_guard(mesh, *cart_opt);
    mpi::world.barrier();
  }

  WHEN("8x8") {
    auto cart_opt = make_cart( {8, 8}, {false, false} );
    if ( cart_opt ) test_merge_guard(mesh, *cart_opt);
    mpi::world.barrier();
  }

  WHEN("16x16") {
    auto cart_opt = make_cart( {16, 16}, {false, false} );
    if ( cart_opt ) test_merge_guard(mesh, *cart_opt);
    mpi::world.barrier();
  }
}

SCENARIO("test merge guard on cartesian topology involving periodic boundary", "[field][mpi]") {
  // TODOL
}


