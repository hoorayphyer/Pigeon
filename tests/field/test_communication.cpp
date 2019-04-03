#include "all_in_one.hpp"
#include "field/communication.cpp"
#include "parallel/mpi++.hpp"

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


// return Ib_bulk and ext
std::tuple<int,int> region_begin_extent ( int region_coord, int mesh_ext, int g ) {
  switch (region_coord) {
  case -1 : return {-g, g};
  case 0 : return {0, mesh_ext - 2 * g};
  case 1 : return {mesh_ext - 2*g, g};
  }
}

// TODO test merge guard
// WHEN("merge guard") {
    //   // only guard cells nonzero values
    //   for ( int i = 0; i < VF::NDim; ++i ) {
    //     double val = field_value( my_coords, i );
    //     for( auto& x : vf[i].data() ) x = val;
    //     for( auto I_bulk : apt::Block<2>({ mesh.bulk_dim(0), mesh.bulk_dim(1) }) ) {
    //       vf[i](I_bulk) = 0.0;
    //     }
    //   }

    //   merge_guard_cells_into_bulk(vf, cart);

    //   apt::Index<2> Ib;
    //   apt::Index<2> ext;
    //   apt::Index<2> region;
    //   // region[ith_dim] can be -1, 0, 1, corresponding to in left guard, in bulk, in right guard
    //   for ( region[1] = -1; region[1] < 2; ++region[1] ) {
    //     std::tie(Ib[1], ext[1]) = region_begin_extent( region[1], mesh.extent()[1], g );
    //     for ( region[0] = -1; region[0] < 2; ++region[0] ) {
    //       std::tie(Ib[0], ext[0]) = region_begin_extent( region[0], mesh.extent()[0], g );
    //       auto neigh = my_coords;
    //       neigh[0] += region[0];
    //       neigh[1] += region[1];

    //       bool at_bdry = false;
    //       for ( int i = 0; i < 2; ++i ) {
    //         if ( 0 == region[i] ) continue;
    //         at_bdry = ( at_bdry || !(cart.shift(i,1)[ 1 == region[i] ]) );
    //       }

    //       // TODO fix these. It has more cases than sync
    //       // for ( int i_fld_dim = 0; i_fld_dim < VF::NDim; ++i_fld_dim ) {
    //       //   double region_val = at_bdry ? 0.0 : field_value( neigh, i_fld_dim );
    //       //   CAPTURE( i_fld_dim, Ib, ext, region, at_bdry );
    //       //   for ( auto I : apt::Block( ext ) ) {
    //       //     REQUIRE( vf[i_fld_dim](I + Ib) == region_val );
    //       //   }
    //       // }
    //     }
    //   }
    // }


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
        at_bdry = ( at_bdry || !(cart.shift(i,1)[ 1 == region[i] ]) );
      }

      for ( int i_fld_dim = 0; i_fld_dim < VF::NDim; ++i_fld_dim ) {
        double region_val = at_bdry ? 0.0 : field_value( neigh, i_fld_dim );
        CAPTURE( i_fld_dim, Ib, ext, region, at_bdry );
        for ( auto I : apt::Block( ext ) ) {
          REQUIRE( vf[i_fld_dim](I + Ib) == region_val );
        }
      }
    }
  }

}

constexpr int g = 1; // guard
constexpr auto mesh = Mesh<2>( {4, 4}, g );

SCENARIO("1x1", "[field][mpi]") {
  auto cart_opt = make_cart( {1, 1}, {false, false} );
  if ( cart_opt ) {
    WHEN("sync guard, which should do nothing in this case") {
      test_sync_guard(mesh, *cart_opt);
    }
  }
  mpi::world.barrier();
}

// TODOL
// SCENARIO("-1x1", "[field][mpi]") {

// }

// TODOL
// SCENARIO("-1x-1", "[field][mpi]") {

// }

SCENARIO("2x1", "[field][mpi]") {
  auto cart_opt = make_cart( {2, 1}, {false, false} );
  if ( cart_opt ) {
    WHEN("sync guard") {
      test_sync_guard(mesh, *cart_opt);
    }
  }
  mpi::world.barrier();
}

SCENARIO("1x2", "[field][mpi]") {
  auto cart_opt = make_cart( {1, 2}, {false, false} );
  if ( cart_opt ) {
    WHEN("sync guard") {
      test_sync_guard(mesh, *cart_opt);
    }
  }
  mpi::world.barrier();

}

SCENARIO("2x2", "[field][mpi]") {
  auto cart_opt = make_cart( {2, 2}, {false, false} );
  if ( cart_opt ) {
    WHEN("sync guard") {
      test_sync_guard(mesh, *cart_opt);
    }
  }
  mpi::world.barrier();
}

SCENARIO("4x4", "[field][mpi]") {
  auto cart_opt = make_cart( {4, 4}, {false, false} );
  if ( cart_opt ) {
    WHEN("sync guard") {
      test_sync_guard(mesh, *cart_opt);
    }
  }
  mpi::world.barrier();
}

SCENARIO("8x8", "[field][mpi]") {
  auto cart_opt = make_cart( {8, 8}, {false, false} );
  if ( cart_opt ) {
    WHEN("sync guard") {
      test_sync_guard(mesh, *cart_opt);
    }
  }
  mpi::world.barrier();
}

SCENARIO("16x16", "[field][mpi]") {
  auto cart_opt = make_cart( {16, 16}, {false, false} );
  if ( cart_opt ) {
    WHEN("sync guard") {
      test_sync_guard(mesh, *cart_opt);
    }
  }
  mpi::world.barrier();
}
