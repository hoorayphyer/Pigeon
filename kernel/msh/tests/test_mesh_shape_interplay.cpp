#include "testfw/testfw.hpp"
#include "msh/mesh_shape_interplay_impl.hpp"
#include "particle/shapef.hpp"

using namespace field;
using namespace msh;

SCENARIO("interpolate component", "[field]") {
  constexpr int DGrid = 2;
  using ShapeF = particle::shapef_t<particle::shape::Cloud_In_Cell>;
  Mesh<DGrid> mesh( {{ {0,16,1}, {0,16,1} }} );
  Field<double,1,2> field(mesh);

  // TODO set up a predefined function, use it to populate the mesh. Then use interpolation and compare with the actual value
}
