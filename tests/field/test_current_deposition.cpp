#include "all_in_one.hpp"
#include "field/communication.cpp"
#include "field/current_deposition.cpp"
#include "kernel/shapef.hpp"
#include "apt/pair.hpp"

using namespace field;
using apt::array;

constexpr int DPtc = 3;

SCENARIO("Testing ShapeRange with 2D grid, deposited field has offsets 0.5, 0.5, 0.5", "[field]") {
  constexpr int DGrid = 2;
  auto test =
    [] ( const auto& shapef, const array<int,DGrid>& cell,
         const array<double,DPtc>& x0_rel, const array<double, DPtc>& dx_rel,
         const array<apt::pair<int>,DGrid>& cells_bounds )
    {
      apt::Vec<double,DPtc> q0_std = {0.0, 0.0, 0.0};
      apt::Vec<double,DPtc> q1_std = {0.0, 0.0, 0.0};

      for ( int i = 0; i < DGrid; ++i ) {
        q0_std[i] = cell[i] + x0_rel[i];
        q1_std[i] = q0_std[i] + dx_rel[i];
      }

      // NOTE need q1 as 1st argument
      const auto[I_b, extent] = impl::deposit_range<DGrid>( q0_std, q1_std, shapef );

      for ( int i = 0; i < DGrid; ++i ) {
        REQUIRE( I_b[i] == cells_bounds[i][0] );
        REQUIRE( extent[i] == cells_bounds[i][1] - cells_bounds[i][0] );
      }

    };


  // GIVEN("Nearest_Grid_Point") {
  //   auto sf = knl::shapef_t<knl::shape::Nearest_Grid_Point>();

  //   WHEN("q0 and q1 affect same X and Y") {
  //     test( sf,
  //           { 50, 50 },
  //           { 0.8, 0.8, 0.0 },
  //           { 0.1, 0.1, 0.0 },
  //           { 50, 51, 50, 51 }
  //           );
  //   }

  //   WHEN("q0 and q1 affect different X but same Y") {
  //     test( sf,
  //           { 50, 50 },
  //           { 0.8, 0.8, 0.0 },
  //           { 0.3, 0.1, 0.0},
  //           { 50, 52, 50, 51 } );
  //   }

  //   WHEN("q0 and q1 affect same X but different Y") {
  //     test( sf,
  //           { 50, 50 },
  //           { 0.8, 0.8, 0.0 },
  //           { 0.1, 0.3, 0.0},
  //           { 50, 51, 50, 52 } );
  //   }

  //   WHEN("q0 and q1 affect same X but different Y") {
  //     test( sf,
  //           { 50, 50 },
  //           { 0.8, 0.8, 0.0 },
  //           { 0.3, 0.3, 0.0},
  //           { 50, 52, 50, 52 } );
  //   }

  // }

  GIVEN("Cloud_In_Cell") {
    auto sf = knl::shapef_t<knl::shape::Cloud_In_Cell>();

    WHEN("q0 and q1 affect same X and Y") {
      test( sf,
            { 50, 50 },
            { 0.3, 0.3, 0.0 },
            { 0.1, 0.1, 0.0},
            { 49, 51, 49, 51 } );
    }
    WHEN("q0 and q1 affect different X but same Y") {
      test( sf,
            { 50, 50 },
            { 0.3, 0.3, 0.0 },
            { 0.4, 0.1, 0.0},
            { 49, 52, 49, 51 } );
    }
    WHEN("q0 and q1 affect same X but different Y") {
      test( sf,
            { 50, 50 },
            { 0.3, 0.3, 0.0 },
            { 0.1, 0.5, 0.0},
            { 49, 51, 49, 52 } );
    }
    WHEN("q0 and q1 affect different X and different Y") {
      test( sf,
            { 50, 50 },
            { 0.3, 0.3, 0.0 },
            { 0.5, 0.5, 0.0},
            { 49, 52, 49, 52 } );

      test( sf,
            { 50, 50 },
            { 0.5, 0.5, 0.0 },
            { 0.5, 0.5, 0.0},
            { 50, 52, 50, 52 } );
    }
  }
}

TEMPLATE_TEST_CASE("Testing dJ in 2D, each cell in bulk has some particles with q0 at the center, and q1 at the top right corner, with q1[2] - q0[2] = 1.0", "[field][mpi]"
                   // NOTE Notation: XxYxZ is the cartesian partition. The cartesian topology is periodic in all directions
                   , (aio::IndexType<-1,-1>)
                   , (aio::IndexType<-2,-1>)
                   , (aio::IndexType<-1,-2>)
                   , (aio::IndexType<-2,-2>)
                   , (aio::IndexType<-4,-4>)
                   // , (aio::IndexType<-8,-8>)
                   ) {
  using ShapeF = knl::shapef_t<knl::shape::Cloud_In_Cell>;
  constexpr ShapeF shapef;
  constexpr int DGrid = 2;
  constexpr int DPtc = 3;

  std::vector<int> cart_dims;
  std::vector<bool> periodic;
  for ( auto i : TestType::get() ) {
    cart_dims.push_back( i > 0 ? i : -i );
    periodic.push_back( true );
  }

  auto cart_opt = aio::make_cart( cart_dims, periodic, mpi::world );

  struct BrutalForce_dJ_Field {
  public:
    using T = double;
  private:
    Field<T,3,DGrid> _data;
    const ShapeF& _shapef;

  public:
    BrutalForce_dJ_Field( const Mesh<DGrid>& mesh, const ShapeF& shapef ) : _data(mesh), _shapef(shapef) {}

    // NOTE q0 and q1 are of type double regardlessly
    // J found by brutal force, which for each ptc, looping over all cells in J and deposit, then immediately integrate to find J
    void deposit( const apt::Vec<double,DPtc>& q0_std, const apt::Vec<double,DPtc>& q1_std ) {
      // simply loop over all cells including guard
      const auto& mesh = _data.mesh();

      // use a "large enough" range to find contributing cells.
      apt::Index<DGrid> Ib;
      apt::Index<DGrid> ext;

      auto int_flr =
        []( T q ) noexcept {
          // Since min(q0_std) = 0.0 by design, min(q_nat) = -1 - r - 0.5, so q_nat + ( support + 3 ) / 2.0 >= 0. We will simply use (support + 1) as the shift
          constexpr auto shift = 1 + ShapeF::support();
          return int(q+shift) - shift; };

      for ( int i_dim = 0; i_dim < DGrid; ++i_dim ) {
        Ib[i_dim] = int_flr( std::min(q0_std[i_dim],q1_std[i_dim]) - 0.5 - ShapeF::support() / 2.0 ) - 1;
        Ib[i_dim] = std::max(Ib[i_dim], -mesh.guard());
        ext[i_dim] = int_flr( std::max(q0_std[i_dim],q1_std[i_dim]) - 0.5 + ShapeF::support() / 2.0 ) + 2;
        ext[i_dim] = std::min(ext[i_dim], mesh.bulk_dim(i_dim) + mesh.guard() );
        ext[i_dim] -= Ib[i_dim];
      }

      apt::array<T,DGrid> s0, s1;
      Field<long double,3,DGrid> W( {ext, 0} ); // use long double to ensure precision
      for( auto I : apt::Block(ext) ) {
        for ( int i = 0; i < DGrid; ++i ) {
          s0[i] = _shapef( I[i] + Ib[i] + 0.5 - q0_std[i] );
          s1[i] = _shapef( I[i] + Ib[i] + 0.5 - q1_std[i] );
        }

        if constexpr ( DGrid == 2 ) {
            W[0](I) = ( s1[0] - s0[0] ) * Wesir( s0[1], s1[1], 1.0, 1.0 );
            W[1](I) = ( s1[1] - s0[1] ) * Wesir( 1.0, 1.0, s0[0], s1[0] );
            W[2](I) = ( q1_std[2] - q0_std[2] ) * Wesir( s0[0], s1[0], s0[1], s1[1] );
          } else {
          W[0](I) = ( s0[0] - s1[0] ) * Wesir( s0[1], s1[1], s0[2], s1[2] );
          W[1](I) = ( s0[1] - s1[1] ) * Wesir( s0[2], s1[2], s0[0], s1[0] );
          W[2](I) = ( s0[2] - s1[2] ) * Wesir( s0[0], s1[0], s0[1], s1[1] );
        }
      }

      // integrate W right away
      for ( int i_dim = 0; i_dim < DGrid; ++i_dim ) {
        for ( auto trI : W.mesh().project(i_dim, W.mesh().origin(), W.mesh().extent() ) ) {
          // J[i+1] = J[i] - W[i], or J[i] = J[i+1] + W[i]. NOTE J is for this particle only
          for ( int n = W.mesh().extent()[i_dim] - 2; n > W.mesh().origin()[i_dim]-1; --n )
            W[i_dim][trI | n ] += W[i_dim][ trI | n+1 ];
        }
      }

      // deposit
      for( auto I : apt::Block(ext) ) {
        for ( int i = 0; i < 3; ++i )
          _data[i](I+Ib) += W[i](I);
      }
    }

    void reduce( int chief, const mpi::Comm& intra ) {
      for ( int i = 0; i < 3; ++i )
        intra.reduce<mpi::IN_PLACE>( mpi::by::SUM, chief, _data[i].data().data(),  _data[i].data().size() );
    }

    auto& integrate( const mpi::CartComm& cart ) {
      merge_guard_cells_into_bulk(_data, cart);
      return _data;
    }

    inline void reset() {
      apt::foreach<0, 3>
        ( []( auto comp ) { //TODOL semantics
            for ( auto& elm : comp.data() ) elm = 0.0;
          }, _data );
    }

  };

  auto test_on_bulk_dims =
    [&shapef] ( const apt::Index<DGrid> bulk_dims, const auto& cart ) {
      THEN("bulk_dims should be large enough") {
        for ( int i = 0; i < 2; ++i )
          REQUIRE( bulk_dims[i] >= 2 * ( (shapef.support() + 3) / 2 ) );
      }

      Standard_dJ_Field<double, 3, DGrid, decltype(shapef)> dJ(bulk_dims, shapef);
      BrutalForce_dJ_Field dJ_bf(dJ.mesh(), shapef);

      dJ.reset();
      dJ_bf.reset();

      const int Nptc_per_cell = 10;

      aio::unif_real<double> unif( -0.9999999, 0.9999999 );

      for ( int j = 0; j < bulk_dims[1]; ++j ) {
        for ( int i = 0; i < bulk_dims[0]; ++i ) {
          for ( int n = 0; n < Nptc_per_cell; ++n ) {
            apt::Vec<double,DPtc> q0( i + std::abs(unif()), j + std::abs(unif()), unif() );
            apt::Vec<double,DPtc> q1( q0[0] + unif(), q0[1] + unif(), unif() );
            dJ.deposit( 1.0, q0, q1 );
            dJ_bf.deposit(q0, q1 );
          }
        }
      }

      auto& Jstd = dJ.integrate( cart );
      auto& Jbf = dJ_bf.integrate( cart );

      const auto& mesh = dJ.mesh();
      const int g = mesh.guard();
      array<bool,DGrid> in_bulk{};
      for ( int comp = 0; comp < 3; ++comp ) {
        for ( int j = -g; j < bulk_dims[1] + g; ++j ) {
          in_bulk[1] = ( j > -1 && j < bulk_dims[1] );
          for ( int i = -g; i < bulk_dims[0] + g; ++i ) {
            in_bulk[0] = ( i > -1 && i < bulk_dims[0] );
            CAPTURE(comp,i,j);
            if ( in_bulk[0] && in_bulk[1] ) {
              REQUIRE( Jstd[comp](apt::Index<2>{i,j}) == Approx(Jbf[comp](apt::Index<2>{i,j})) );
            } else {
              // in guard cells it should be already cleared out
              REQUIRE( Jstd[comp](apt::Index<2>{i,j}) == 0.0 );
              REQUIRE( Jbf[comp](apt::Index<2>{i,j}) == 0.0 );
            }
          }
        }
      }
    };

  if ( cart_opt ) {
    test_on_bulk_dims({4,4}, *cart_opt);
    test_on_bulk_dims({8,8}, *cart_opt);
  }

  mpi::world.barrier();
}
