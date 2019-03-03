#include "particle/depositer.cpp"
#include "kernel/grid.hpp"
#include "kernel/shapef.hpp"
#include "catch2/catch.hpp"

using namespace particle;
using namespace esirkepov;

constexpr int DPtc = 3;

SCENARIO("Testing ShapeRange with 2D grid", "[particle]") {
  constexpr int DGrid = 2;
  using Grid = knl::Grid<DGrid, double>;

  auto test
    = []( Grid grid, auto shapef, std::array<int,DGrid> cell,
          std::array<double,DPtc> x0_rel,
          std::array<double, DPtc> dx_rel,
          std::array<std::array<int,2>,DGrid> cells_bounds ) {
        for( int i = 0; i < DGrid; ++i )
          cells_bounds[i][1] -= cells_bounds[i][0]; // get the stride

        apt::Vec<double,DPtc> q1_abs = {0.0, 0.0, 0.0};
        std::array<double,DPtc> dq_abs = {0.0, 0.0, 0.0};

        for ( int i = 0; i < DGrid; ++i ) {
          dq_abs[i] = dx_rel[i] * grid[i].delta;
          q1_abs[i] = grid[i].absc(cell[i], x0_rel[i]) + dq_abs[i];
        }
        CAPTURE(q1_abs);



        auto sr = make_shape_range(q1_abs, apt::Vec<double,DPtc>{dq_abs}, grid, shapef ); // NOTE need q1 as 1st argument

        REQUIRE( sr.end() - sr.begin() == cells_bounds[0][1] * cells_bounds[1][1] );
        int I = 0;
        std::array<int, DGrid> ijk{};
        for (auto[I_n,W_n] : sr) {
          ijk = { I % cells_bounds[0][1] + cells_bounds[0][0],
                  I / cells_bounds[0][1] + cells_bounds[1][0]
          };
          REQUIRE( I_n == ijk );
          ++I;
        }
      };

  GIVEN("Nearest_Grid_Point") {
    auto sf = knl::shapef_t<knl::shape::Nearest_Grid_Point>();
    Grid grid{ {{
                 { 0.0, 1.0, 0, 100 },
                 { 0.0, 1.0, 0, 100 }
                   }} };


    WHEN("q0 and q1 affect same X and Y") {
      test( grid, sf,
            { 50, 50 },
            { 0.8, 0.8, 0.0 },
            { 0.1, 0.1, 0.0},
            { 50, 51, 50, 51 } );
    }

    WHEN("q0 and q1 affect different X but same Y") {
      test( grid, sf,
            { 50, 50 },
            { 0.8, 0.8, 0.0 },
            { 0.3, 0.1, 0.0},
            { 50, 52, 50, 51 } );
    }

    WHEN("q0 and q1 affect same X but different Y") {
      test( grid, sf,
            { 50, 50 },
            { 0.8, 0.8, 0.0 },
            { 0.1, 0.3, 0.0},
            { 50, 51, 50, 52 } );
    }

    WHEN("q0 and q1 affect same X but different Y") {
      test( grid, sf,
            { 50, 50 },
            { 0.8, 0.8, 0.0 },
            { 0.3, 0.3, 0.0},
            { 50, 52, 50, 52 } );
    }


  }

  GIVEN("Cloud_In_Cell") {
    auto sf = knl::shapef_t<knl::shape::Cloud_In_Cell>();
    Grid grid{ {{
                 { 1.0, 2.0, 1, 100 },
                 { 1.0, 2.0, 1, 100 }
                   }} };

    WHEN("q0 and q1 affect same X and Y") {
      test( grid, sf,
            { 50, 50 },
            { 0.3, 0.3, 0.0 },
            { 0.1, 0.1, 0.0},
            { 49, 51, 49, 51 } );
    }

    WHEN("q0 and q1 affect different X but same Y") {
      test( grid, sf,
            { 50, 50 },
            { 0.3, 0.3, 0.0 },
            { 0.4, 0.1, 0.0},
            { 49, 52, 49, 51 } );
    }
    WHEN("q0 and q1 affect same X but different Y") {
      test( grid, sf,
            { 50, 50 },
            { 0.3, 0.3, 0.0 },
            { 0.1, 0.5, 0.0},
            { 49, 51, 49, 52 } );
    }
    WHEN("q0 and q1 affect different X and different Y") {
      test( grid, sf,
            { 50, 50 },
            { 0.3, 0.3, 0.0 },
            { 0.5, 0.5, 0.0},
            { 49, 52, 49, 52 } );
    }
  }


}
