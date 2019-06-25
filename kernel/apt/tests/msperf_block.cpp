#include "testfw/testfw.hpp"
#include "apt/block.hpp"
#include "timer/timer.hpp"

using namespace apt;

SCENARIO("Measure performance of Block", "[apt]") {

  constexpr Index<3> extent {128, 128, 128};
  Block<3> block (extent);

  const int N = 1000;
  tmr::Timestamp t1;
  for ( int dummy = 0; dummy < N; ++dummy ) {

    for ( const auto& I : block ) {
      auto x = I[0] + I[1] * I[2];
    }

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
  std::cout << dur.val() << std::endl;
  WARN("Timing Block of " << extent << " = " << dur.val() / N << " " << dur.unit() );

}
