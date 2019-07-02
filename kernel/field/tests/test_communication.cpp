#include "testfw/testfw.hpp"
#include "field/communication_impl.hpp"
#include "mpipp/mpi++.hpp"
#include <algorithm> // for std::min
#include "debug/nan.hpp"

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

void test_sync_guard ( const Mesh<2>& mesh, const mpi::CartComm& cart ) {
  const auto[ my_coords, dims, periodic] = cart.coords_dims_periodic();
  VF vf( mesh );

  // only bulk has nonzero values
  for ( int i = 0; i < VF::NDim; ++i ) {
    double val = field_value( my_coords, i );
    for( auto I_bulk : apt::Block<2>({ mesh.bulk_dim(0), mesh.bulk_dim(1) }) ) {
      vf[i](I_bulk) = val;
    }
  }

  REQUIRE( debug::num_nan(vf) == 0 );

  sync_guard_cells_from_bulk(vf, cart);

  REQUIRE( debug::num_nan(vf) == 0 );

  // return Ib_bulk and ext
  auto region_begin_extent =
    []( int lcr, int mesh_ext, int g ) -> std::tuple<int,int>
    {
     switch (lcr) {
     case -1 : return {-g, g};
     case 0 : return {0, mesh_ext - 2 * g};
     case 1 : return {mesh_ext - 2*g, g};
     }
    };

  apt::Index<2> Ib;
  apt::Index<2> ext;
  apt::Index<2> lcr;
  // lcr[ith_dim] can be -1, 0, 1, corresponding to in left guard, in bulk, in right guard
  for ( lcr[1] = -1; lcr[1] < 2; ++lcr[1] ) {
    std::tie(Ib[1], ext[1]) = region_begin_extent( lcr[1], mesh.extent()[1], mesh.guard() );
    for ( lcr[0] = -1; lcr[0] < 2; ++lcr[0] ) {
      std::tie(Ib[0], ext[0]) = region_begin_extent( lcr[0], mesh.extent()[0], mesh.guard() );
      auto neigh = my_coords;
      for ( int i_dim = 0; i_dim < 2; ++i_dim )
        neigh[i_dim] = ( my_coords[i_dim] + lcr[i_dim] + dims[i_dim] ) % dims[i_dim];

      bool at_bdry = false;
      for ( int i_dim = 0; i_dim < 2; ++i_dim ) {
        if ( 0 == lcr[i_dim] ) continue;
        at_bdry = ( at_bdry || !(cart.shift(i_dim)[ 1 == lcr[i_dim] ]) );
      }

      INFO("checking bulk areas in sync_guard");
      for ( int i_fld_dim = 0; i_fld_dim < VF::NDim; ++i_fld_dim ) {
        double region_val = at_bdry ? 0.0 : field_value( neigh, i_fld_dim );
        CAPTURE( mpi::world.rank(), i_fld_dim, Ib, ext, lcr, at_bdry );
        for ( auto I : apt::Block( ext ) )
          REQUIRE( vf[i_fld_dim](I + Ib) == region_val );
      }
    }
  }
}

void test_merge_guard ( const Mesh<2>& mesh, const mpi::CartComm& cart ) {
  const auto[ my_coords, dims, periodic] = cart.coords_dims_periodic();
  VF vf( mesh );

  // all cells have nonzero values
  for ( int i = 0; i < VF::NDim; ++i ) {
    double val = field_value( my_coords, i );
    for( auto& x : vf[i].data() ) x = val;
  }

  REQUIRE( debug::num_nan(vf) == 0 );

  merge_guard_cells_into_bulk(vf, cart);

  REQUIRE( debug::num_nan(vf) == 0 );
  { // Check bulk value

    // return Ib_bulk and ext
    auto region_begin_extent =
      []( int lcr, int mesh_ext, int g ) -> std::tuple<int,int>
      {
       switch (lcr) {
       case -1 : return {0, g};
       case 0 : return {g, mesh_ext - 4 * g};
       case 1 : return {mesh_ext - 3*g, g};
       }
      };

    apt::Index<2> Ib;
    apt::Index<2> ext;
    apt::Index<2> lcr;
    // lcr[ith_dim] can be -1, 0, 1, corresponding to left part in bulk of width guard, bulk, right part in bulk of width guard
    for ( lcr[1] = -1; lcr[1] < 2; ++lcr[1] ) {
      std::tie(Ib[1], ext[1]) = region_begin_extent( lcr[1], mesh.extent()[1], mesh.guard() );
      for ( lcr[0] = -1; lcr[0] < 2; ++lcr[0] ) {
        std::tie(Ib[0], ext[0]) = region_begin_extent( lcr[0], mesh.extent()[0], mesh.guard() );

        for ( int i_fld_dim = 0; i_fld_dim < VF::NDim; ++i_fld_dim ) {
          double region_val = 0.0;
          // RATIONALE for (regx, regy), contributing neighbors are set of { ( regx, 0 ), ( regx, regy ), (0, regy), (0,0) } with duplicates all dropped. Note (0,0) is self. When reg_i = 1 or -1 corresponds to topology boundary, that dimension is ignored.
          apt::Index<2> contrib;
          for ( contrib[1] = std::min(lcr[1],0); contrib[1] <= std::max(lcr[1],0); ++contrib[1] ) {
            if ( contrib[1] != 0 && !cart.shift(1)[1 == contrib[1]] ) continue; // check if at cartesian topology end
            for ( contrib[0] = std::min(lcr[0],0); contrib[0] <= std::max(lcr[0],0); ++contrib[0] ) {
              if ( contrib[0] != 0 && !cart.shift(0)[1 == contrib[0]] ) continue; // check if at cartesian topology end
              auto neigh = my_coords;
              for ( int i_dim = 0; i_dim < 2; ++i_dim )
                neigh[i_dim] = ( my_coords[i_dim] + contrib[i_dim] + dims[i_dim] ) % dims[i_dim];
              region_val += field_value( neigh, i_fld_dim );
            }
          }

          INFO("checking bulk areas in merge_guard");
          CAPTURE( mpi::world.rank(), i_fld_dim, Ib, ext, lcr  );
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


TEMPLATE_TEST_CASE("field commuication on cartesian topology", "[field][mpi]"
                   // nonperiodic
                   , (aio::IndexType<1,1>)
                   , (aio::IndexType<2,1>)
                   , (aio::IndexType<1,2>)
                   , (aio::IndexType<2,2>)
                   , (aio::IndexType<4,4>)
                   , (aio::IndexType<8,8>)
                   , (aio::IndexType<16,16>)
                   // periodic
                   , (aio::IndexType<-1,1>)
                   , (aio::IndexType<1,-1>)
                   , (aio::IndexType<-1,-1>)
                   , (aio::IndexType<-2,1>)
                   , (aio::IndexType<2,-1>)
                   , (aio::IndexType<-2,-2>)
                   ) {
  constexpr auto CartTopo = TestType::get();
  CAPTURE(CartTopo);
  constexpr auto mesh = Mesh<2>( {4, 4}, 1 );
  THEN("left and right deposited areas should not overlap") {
    for ( int i = 0; i < 2; ++i )
      REQUIRE( mesh.extent()[i] >= 4 * mesh.guard() );
  }

  std::vector<int> cart_dims;
  std::vector<bool> periodic;
  for ( auto i : CartTopo ) {
    bool is_neg = (i < 0);
    cart_dims.push_back( is_neg ? -i : i );
    periodic.push_back( is_neg );
  }

  // NOTE somehow using same cart_opt across sync and merge tests cause the program to crash
  WHEN("syncing guard") {
    auto cart_opt = aio::make_cart( cart_dims, periodic, mpi::world );
    if ( cart_opt ) test_sync_guard(mesh, *cart_opt);
  }
  mpi::world.barrier();

  WHEN("merging guard") {
    auto cart_opt = aio::make_cart( cart_dims, periodic, mpi::world );
    if ( cart_opt) test_merge_guard(mesh, *cart_opt);
  }
  mpi::world.barrier();

}
