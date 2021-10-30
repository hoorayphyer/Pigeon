#include "msh/mesh_shape_interplay_impl.hpp"
#include "particle/shapef.hpp"
#include "testfw/testfw.hpp"
#include "timer/timer.hpp"

using namespace msh;

SCENARIO("Measure performance of interpolate", "[msh]") {
  using ShapeF = particle::shapef_t<particle::shape::Cloud_In_Cell>;

  constexpr int guard = 1;  // TODO this should depend on ShapeF
  const field::Mesh<2> mesh(apt::make_range<2>({}, {128, 128}, guard));

  field::Field<double, 3, 2> E(mesh);
  E.set_offset(0, {INSITU, MIDWAY});
  E.set_offset(1, {MIDWAY, INSITU});
  E.set_offset(2, {MIDWAY, MIDWAY});

  aio::unif_real<double> unif;
  for (int C = 0; C < 3; ++C) {
    for (auto& x : E[C].data()) x = 1000 * (2 * unif() - 1);
  }

  auto shapef = ShapeF();

  const int N = 1;
  std::vector<double> qs(2 * N);  // 2 = DGrid
  for (auto& x : qs) x = unif() * 127.9;

  apt::array<double, 2> q_std{50, 50};

  tmr::Timestamp t1;

  for (int n = 0; n < N; ++n) {
    q_std = {qs[2 * n], qs[2 * n + 1]};

    auto E_itpl = msh::interpolate(E, q_std, shapef);
    // for ( int k = 0; k < extent[2]; ++k ) {
    //   for ( int j = 0; j < extent[1]; ++j ) {
    //     for ( int i = 0; i < extent[0]; ++i ) {
    //       auto x = i + j * k;
    //     }
    //   }
    // }
  }

  auto dur = t1.lapse();
  dur.in_units_of("ns");
  WARN("Timing Interpolating fields of " << apt::range::size(mesh.range())
                                         << " = " << dur.val() / N << " "
                                         << dur.unit());
}
