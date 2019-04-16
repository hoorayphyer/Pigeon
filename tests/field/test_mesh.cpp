#include "all_in_one.hpp"
#include "field/mesh.hpp"

using namespace field;

SCENARIO("Test Mesh", "[field][mesh]") {
  constexpr int D = 2;
  apt::Index<D> bulk_ext { 16, 256 };
  Mesh mesh( bulk_ext, 2 );

  REQUIRE(mesh.guard() == 2);
  REQUIRE(mesh.origin() == apt::Index<D>{-2, -2});
  REQUIRE(mesh.extent() == apt::Index<D>{20, 260});
  for ( int i = 0; i < D; ++i ) {
    REQUIRE(mesh.bulk_dim(i) == bulk_ext[i] );
  }

  REQUIRE( mesh.stride(0) == 1 );
  REQUIRE( mesh.stride(1) == 20 );
}

SCENARIO("Test LinearizedIndex", "[field][mesh]") {
  constexpr int D = 2;
  apt::Index<D> bulk_ext { 128, 256 };
  Mesh mesh( bulk_ext, 1 );

  apt::Index<D> stride { 1, 130 };
  aio::unif_int<unsigned int> dist(0, 255);

  int N = 100*1000;
  while( N-- ) {
    apt::Index<D> data { dist(), dist() }; // NOTE the index can go out of bounds, but it is fine.
    auto index = mesh.linearized_index_of_whole_mesh(data);
    REQUIRE( index == (data[0] + 1) * stride[0] +  (data[1] + 1) * stride[1] );
  }

}

SCENARIO("Test ProjBlock and TransIndex", "[field][mesh]") {
  constexpr int D = 2;
  apt::Index<D> bulk_ext { 128, 256 };
  Mesh mesh( bulk_ext, 1 );

  apt::Index<D> stride { 1, 130 };
  aio::unif_int<unsigned int> dist(0, 255);

  int N = 100*1000;
  std::cout << "Caution: big loop is going on. Patience please!" << std::endl;
  while( N-- ) {
    apt::Index<D> I_bulk_begin { dist(), dist() }; // NOTE the index can go out of bounds, but it is fine.
    apt::Index<D> extent { dist(), dist() };
    int i_dim = dist() % D;
    auto pb = mesh.project(i_dim, I_bulk_begin, extent);
    auto itr = pb.begin();
    for ( int tr = I_bulk_begin[1-i_dim]; tr < I_bulk_begin[1-i_dim] + extent[1-i_dim]; ++tr ) {
      auto trI = *(itr++);
      for ( int n = I_bulk_begin[i_dim]; n < I_bulk_begin[i_dim] + extent[i_dim]; ++n ) {
        int x, y;
        if ( 0 == i_dim ) {
          x = n;
          y = tr;
        } else {
          x = tr;
          y = n;
        }
        REQUIRE(  (trI | n) == (x + 1) * stride[0] +  (y+1) * stride[1] );
      }
    }
  }

}

// SCENARIO("Test mesh squeeze", "[field][mesh]") {
//   constexpr int D = 2;
//   apt::Index<D> bulk_ext { 128, 256 };
//   Mesh mesh( bulk_ext, 1 );
//   auto mesh_sq0 = mesh.squeeze(0);
//   REQUIRE( mesh_sq0.extent()[0] == 1 + 2 * mesh.guard() );
//   REQUIRE( mesh_sq0.extent()[1] == mesh.extent()[1] );

//   auto mesh_sq1 = mesh.squeeze(1);
//   REQUIRE( mesh_sq1.extent()[1] == 1 + 2 * mesh.guard() );
//   REQUIRE( mesh_sq1.extent()[0] == mesh.extent()[0] );
// }
