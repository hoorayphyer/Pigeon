#include "apt/pair.hpp"
#include "msh/current_deposition_impl.hpp"
#include "particle/shapef.hpp"
#include "silopp/silo++.hpp"
#include "testfw/testfw.hpp"

using namespace field;
using namespace msh;
using apt::array;

constexpr int DPtc = 3;

SCENARIO(
    "Testing ShapeRange with 2D grid, deposited field has offsets 0.5, 0.5, "
    "0.5",
    "[field]") {
  constexpr int DGrid = 2;
  auto test = [](const auto& shapef, const array<int, DGrid>& cell,
                 const array<double, DPtc>& x0_rel,
                 const array<double, DPtc>& dx_rel,
                 const array<apt::pair<int>, DGrid>& cells_bounds) {
    apt::array<double, DPtc> q0_std = {0.0, 0.0, 0.0};
    apt::array<double, DPtc> q1_std = {0.0, 0.0, 0.0};

    for (int i = 0; i < DGrid; ++i) {
      q0_std[i] = cell[i] + x0_rel[i];
      q1_std[i] = q0_std[i] + dx_rel[i];
    }

    // NOTE need q1 as 1st argument
    const auto [I_b, extent] =
        impl::deposit_range<DGrid>(q0_std, q1_std, shapef);

    for (int i = 0; i < DGrid; ++i) {
      REQUIRE(I_b[i] == cells_bounds[i][0]);
      REQUIRE(extent[i] == cells_bounds[i][1] - cells_bounds[i][0]);
    }
  };

  // GIVEN("Nearest_Grid_Point") {
  //   auto sf = particle::shapef_t<particle::shape::Nearest_Grid_Point>();

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
    auto sf = particle::shapef_t<particle::shape::Cloud_In_Cell>();

    WHEN("q0 and q1 affect same X and Y") {
      test(sf, {50, 50}, {0.3, 0.3, 0.0}, {0.1, 0.1, 0.0}, {49, 51, 49, 51});
    }
    WHEN("q0 and q1 affect different X but same Y") {
      test(sf, {50, 50}, {0.3, 0.3, 0.0}, {0.4, 0.1, 0.0}, {49, 52, 49, 51});
    }
    WHEN("q0 and q1 affect same X but different Y") {
      test(sf, {50, 50}, {0.3, 0.3, 0.0}, {0.1, 0.5, 0.0}, {49, 51, 49, 52});
    }
    WHEN("q0 and q1 affect different X and different Y") {
      test(sf, {50, 50}, {0.3, 0.3, 0.0}, {0.5, 0.5, 0.0}, {49, 52, 49, 52});

      test(sf, {50, 50}, {0.5, 0.5, 0.0}, {0.5, 0.5, 0.0}, {50, 52, 50, 52});
    }
  }
}

// NOTE q0 and q1 are of type double regardlessly
// J found by brutal force, which for each ptc, looping over all cells in J and
// deposit, then immediately integrate to find J
template <typename RJ, int DGrid, int DPtc, typename ShapeF>
void direct_deposit(Field<RJ, 3, DGrid>& J, const ShapeF& shapef,
                    const apt::array<double, DPtc>& q0_std,
                    const apt::array<double, DPtc>& q1_std) {
  // simply loop over all cells including guard
  const auto& mesh = J.mesh();
  assert(apt::is_margin_uniform(mesh.range()));
  const int guard = mesh.range()[0].margin()[LFT];

  // use a "large enough" range to find contributing cells.
  apt::Index<DGrid> Ib;
  apt::Index<DGrid> ext;

  auto int_flr = [](auto q) noexcept {
    // Since min(q0_std) = 0.0 by design, min(q_nat) = -1 - r - 0.5, so q_nat +
    // ( support + 3 ) / 2.0 >= 0. We will simply use (support + 1) as the shift
    constexpr auto shift = 1 + ShapeF::support();
    return int(q + shift) - shift;
  };

  for (int i_dim = 0; i_dim < DGrid; ++i_dim) {
    Ib[i_dim] = int_flr(std::min(q0_std[i_dim], q1_std[i_dim]) - 0.5 -
                        ShapeF::support() / 2.0) -
                1;
    Ib[i_dim] = std::max(Ib[i_dim], -guard);
    ext[i_dim] = int_flr(std::max(q0_std[i_dim], q1_std[i_dim]) - 0.5 +
                         ShapeF::support() / 2.0) +
                 2;
    ext[i_dim] = std::min(ext[i_dim], mesh.range(i_dim).size() + guard);
    ext[i_dim] -= Ib[i_dim];
  }

  apt::array<double, DGrid> s0, s1;
  apt::array<apt::Range, DGrid> range;
  for (int i = 0; i < DGrid; ++i) range[i] = {0, ext[i], 0};
  Field<RJ, 3, DGrid> W{range};
  for (auto I : apt::Block({}, ext)) {  // POLEDANCE use Ib and Ie for W
    for (int i = 0; i < DGrid; ++i) {
      s0[i] = shapef(I[i] + Ib[i] + 0.5 - q0_std[i]);
      s1[i] = shapef(I[i] + Ib[i] + 0.5 - q1_std[i]);
    }

    if constexpr (DGrid == 2) {
      W[0](I) = (s1[0] - s0[0]) * Wesir(s0[1], s1[1], 1.0, 1.0);
      W[1](I) = (s1[1] - s0[1]) * Wesir(1.0, 1.0, s0[0], s1[0]);
      W[2](I) = (q1_std[2] - q0_std[2]) * Wesir(s0[0], s1[0], s0[1], s1[1]);
    } else {
      W[0](I) = (s0[0] - s1[0]) * Wesir(s0[1], s1[1], s0[2], s1[2]);
      W[1](I) = (s0[1] - s1[1]) * Wesir(s0[2], s1[2], s0[0], s1[0]);
      W[2](I) = (s0[2] - s1[2]) * Wesir(s0[0], s1[0], s0[1], s1[1]);
    }
  }

  {  // integrate W right away. J[i+1] = J[i] - W[i], or J[i] = J[i+1] + W[i],
     // where J is for this particle only, and we will store it in W

    const auto& r = W.mesh().range();
    for (int i_dim = 0; i_dim < DGrid; ++i_dim) {
      for (const auto& trI : apt::project_out(i_dim, apt::range::far_begin(r),
                                              apt::range::far_end(r))) {
        for (apt::Longidx n(i_dim, r[i_dim].full_size() - 2);
             n > r[i_dim].far_begin() - 1; --n)
          W[i_dim](trI + n) += W[i_dim](trI + (n + 1));
      }
    }

    // NOTE the following alternative implementation uses traditional loops
    // if constexpr ( DGrid == 2 ) {
    //     for ( int j = 0; j < ext[1]; ++j ) {
    //       for ( int i = ext[0] - 2; i > -1; --i )
    //         W[0]({i,j}) += W[0]( {i+1, j} );
    //     }

    //     for ( int i = 0; i < ext[0]; ++i ) {
    //       for ( int j = ext[1] - 2; j > -1; --j )
    //         W[1]({i,j}) += W[1]( {i, j+1} );
    //     }
    //   }
  }

  // deposit
  // POLEDANCE
  for (auto I : apt::Block({}, ext)) {
    for (int i = 0; i < 3; ++i) J[i](I + Ib) += W[i](I);
  }
}

TEMPLATE_TEST_CASE(
    "Testing deposition in 2D against alternative implementation "
    "BrutalForce_dJ_Field",
    "[field]"
    // , (aio::IndexType<4,4>)
    // , (aio::IndexType<8,8>)
    // , (aio::IndexType<16,16>)
    // , (aio::IndexType<32,32>)
    ,
    (aio::IndexType<64, 64>)
    // , (aio::IndexType<128,128>)
    // , (aio::IndexType<256,256>)
    // , (aio::IndexType<512,512>)
) {
  using ShapeF = particle::shapef_t<particle::shape::Cloud_In_Cell>;
  constexpr ShapeF shapef;
  constexpr int DGrid = 2;
  constexpr int DPtc = 3;

  apt::Index<DGrid> bulk_dims{};
  for (int i = 0; i < DGrid; ++i) bulk_dims[i] = TestType::get()[i];

  THEN("bulk_dims should be large enough") {
    for (int i = 0; i < DGrid; ++i)
      REQUIRE(bulk_dims[i] >= 2 * ((shapef.support() + 3) / 2));
  }

  using RJ = double;
  constexpr auto guard = (ShapeF::support() + 3) / 2;

  apt::array<apt::Range, DGrid> range;
  for (int i = 0; i < DGrid; ++i) range[i] = {0, bulk_dims[i], guard};

  Field<RJ, 3, DGrid> J_std{range};
  Field<RJ, 3, DGrid> J_direct{range};

  J_std.reset();
  J_direct.reset();

  const int Nptc_per_cell = 10;

  aio::unif_real<double> unif(-0.9999999, 0.9999999);

  for (int j = 0; j < bulk_dims[1]; ++j) {
    for (int i = 0; i < bulk_dims[0]; ++i) {
      for (int n = 0; n < Nptc_per_cell; ++n) {
        apt::array<double, DPtc> q0{i + std::abs(unif()), j + std::abs(unif()),
                                    unif()};
        apt::array<double, DPtc> q1{q0[0] + unif(), q0[1] + unif(), unif()};
        deposit(J_std, 1.0, shapef, q0, q1);
        direct_deposit(J_direct, shapef, q0, q1);
      }
    }
  }

  // integrate(J_std);

  const auto mesh_size = J_std.mesh().stride().back();

  for (int iJ = 0; iJ < 3; ++iJ) {
    const auto& j_std = J_std[iJ].data();
    const auto& j_direct = J_direct[iJ].data();
    for (int i = 0; i < mesh_size; ++i) {
      REQUIRE(j_std[i] == Approx(j_direct[i]).margin(1e-12));
    }
  }

  if (false) {  // write to silo to visualize the result
    using namespace silo;
    auto dbfile = open("test_current_deposition.silo", Mode::Write);
    {                               // put a quad mesh
      apt::array<int, DGrid> dims;  // number of silo nodes
      std::vector<int> lo_ofs(DGrid);
      std::vector<int> hi_ofs(DGrid);
      for (int i = 0; i < DGrid; ++i) {
        dims[i] = bulk_dims[0] + 2 * guard + 1;
        lo_ofs[i] = guard;
        hi_ofs[i] = guard;
      }

      std::vector<std::vector<double>> coords(dims.size());
      for (int i = 0; i < 2; ++i) {
        coords[i].resize(dims[i]);
        for (int j = 0; j < dims[i]; ++j) {
          coords[i][j] = j - lo_ofs[i];
        }
      }

      OptList optlist;
      // optlist[Opt::LO_OFFSET] = lo_ofs;
      // optlist[Opt::HI_OFFSET] = hi_ofs;

      dbfile.put_mesh("mesh", coords, MeshType::Rect, optlist);
    }
    {  // put J fields
      std::vector<int> data_dims(DGrid);
      for (int i = 0; i < DGrid; ++i)
        data_dims[i] = J_std.mesh().range()[i].full_size();

      dbfile.put_var("J1_std", "mesh", J_std[0].data().data(), data_dims);
      dbfile.put_var("J1_direct", "mesh", J_direct[0].data().data(), data_dims);
      dbfile.put_var("J2_std", "mesh", J_std[1].data().data(), data_dims);
      dbfile.put_var("J2_direct", "mesh", J_direct[1].data().data(), data_dims);
    }
  }
}
