#pragma once

#include "apt/grid.hpp"
#include "field/field.hpp"

namespace mpi {
struct CartComm;
}

namespace field {
template <typename R, int DGrid, typename RJ>
struct LogSphericalSolver {
  // NOTE Assume second order difference scheme
  // NOTE Ignore setting of range for now
 private:
  R _fourpi{};
  R _surface{};  // the cell containing _surface is the last cell below surface.
                 // Using R here so that each local grid knows if it's at
                 // surface boundary
  R _outer{};
  int _op_inv_precision = 4;
  R _alpha = 1.0;

 public:
  LogSphericalSolver(R four_pi, R alpha, int op_inv_precision, R surface,
                     R outer)
      : _fourpi(four_pi),
        _surface(surface),
        _outer(outer),
        _op_inv_precision(op_inv_precision),
        _alpha(alpha) {}

  void operator()(Field<R, 3, DGrid>& E, Field<R, 3, DGrid>& B,
                  Field<RJ, 3, DGrid>& J, const apt::Grid<R, DGrid>& grid,
                  const mpi::CartComm& cart, int timestep, R dt) const;

  static constexpr int min_guard(int iteration_in_inverting_operator) {
    return 2 + iteration_in_inverting_operator;
  }
};

}  // namespace field
