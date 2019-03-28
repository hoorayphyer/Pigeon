#include "field/current_deposition.cpp"
#include "kernel/shapef.hpp"
#include "apt/print.hpp"
#include "apt/pair.hpp"
#include "catch2/catch.hpp"
#include <iostream>

using namespace field;
using apt::array;

constexpr int DPtc = 3;

template < int DGrid, typename ShapeF >
void test ( const ShapeF& shapef,
            const array<int,DGrid>& cell,
            const array<double,DPtc>& x0_rel,
            const array<double, DPtc>& dx_rel,
            const array<apt::pair<int>,DGrid>& cells_bounds ) {
  apt::Vec<double,DPtc> q0_std = {0.0, 0.0, 0.0};
  apt::Vec<double,DPtc> q1_std = {0.0, 0.0, 0.0};

  for ( int i = 0; i < DGrid; ++i ) {
    q0_std[i] = cell[i] + x0_rel[i];
    q1_std[i] = q0_std[i] + dx_rel[i];
  }

  // NOTE need q1 as 1st argument
  const auto[I_b, extent, sep0_b, sep1_b] = impl::set_up_for_depositing_dJ<DGrid>( q0_std, q1_std, shapef );

  for ( int i = 0; i < DGrid; ++i ) {
    REQUIRE( I_b[i] == cells_bounds[i][0] );
    REQUIRE( extent[i] == cells_bounds[i][1] - cells_bounds[i][0] );
    REQUIRE( sep0_b[i] == Approx(cells_bounds[i][0] + 0.5 - cell[i] - x0_rel[i]) );
    REQUIRE( sep1_b[i] - sep0_b[i] == Approx(-dx_rel[i]) );
  }

};

SCENARIO("Testing ShapeRange with 2D grid, deposited field has offsets 0.5, 0.5, 0.5", "[particle]") {
  constexpr int DGrid = 2;

  GIVEN("Nearest_Grid_Point") {
    auto sf = knl::shapef_t<knl::shape::Nearest_Grid_Point>();

    WHEN("q0 and q1 affect same X and Y") {
      test<DGrid>( sf,
                   { 50, 50 },
                   { 0.8, 0.8, 0.0 },
                   { 0.1, 0.1, 0.0 },
                   { 50, 51, 50, 51 }
                   );
    }

    WHEN("q0 and q1 affect different X but same Y") {
      test<DGrid>( sf,
                   { 50, 50 },
                   { 0.8, 0.8, 0.0 },
                   { 0.3, 0.1, 0.0},
                   { 50, 52, 50, 51 } );
    }

    WHEN("q0 and q1 affect same X but different Y") {
      test<DGrid>( sf,
                   { 50, 50 },
                   { 0.8, 0.8, 0.0 },
                   { 0.1, 0.3, 0.0},
                   { 50, 51, 50, 52 } );
    }

    WHEN("q0 and q1 affect same X but different Y") {
      test<DGrid>( sf,
                   { 50, 50 },
                   { 0.8, 0.8, 0.0 },
                   { 0.3, 0.3, 0.0},
                   { 50, 52, 50, 52 } );
    }

  }

  GIVEN("Cloud_In_Cell") {
    auto sf = knl::shapef_t<knl::shape::Cloud_In_Cell>();

    WHEN("q0 and q1 affect same X and Y") {
      test<DGrid>( sf,
                   { 50, 50 },
                   { 0.3, 0.3, 0.0 },
                   { 0.1, 0.1, 0.0},
                   { 49, 51, 49, 51 } );
    }

    WHEN("q0 and q1 affect different X but same Y") {
      test<DGrid>( sf,
                   { 50, 50 },
                   { 0.3, 0.3, 0.0 },
                   { 0.4, 0.1, 0.0},
                   { 49, 52, 49, 51 } );
    }
    WHEN("q0 and q1 affect same X but different Y") {
      test<DGrid>( sf,
                   { 50, 50 },
                   { 0.3, 0.3, 0.0 },
                   { 0.1, 0.5, 0.0},
                   { 49, 51, 49, 52 } );
    }
    WHEN("q0 and q1 affect different X and different Y") {
      test<DGrid>( sf,
                   { 50, 50 },
                   { 0.3, 0.3, 0.0 },
                   { 0.5, 0.5, 0.0},
                   { 49, 52, 49, 52 } );
    }
  }

}
