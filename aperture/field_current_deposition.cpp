#include "field/current_deposition.cpp"
#include "traits.hpp"
namespace field {
  using namespace traits;
  template class Standard_dJ_Field< real_dj_t, 3, DGrid, ShapeF >;
}
