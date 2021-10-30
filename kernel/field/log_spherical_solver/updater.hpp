#ifndef _FIELD_LOG_SPHERICAL_UPDATER_HPP_
#define _FIELD_LOG_SPHERICAL_UPDATER_HPP_

#include "field/action.hpp"

namespace field {
template <typename R, int DGrid, typename RJ>
struct LogSphericalSolver : public Action<R, DGrid, RJ> {
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
  constexpr auto& set_fourpi(R x) noexcept {
    _fourpi = x;
    return *this;
  }
  constexpr auto& set_op_inv_precision(int n) noexcept {
    _op_inv_precision = n;
    return *this;
  }
  constexpr auto& set_alpha(R a) noexcept {
    _alpha = a;
    return *this;
  }
  constexpr auto& set_surface(R a) noexcept {
    _surface = a;
    return *this;
  }
  constexpr auto& set_outer(R a) noexcept {
    _outer = a;
    return *this;
  }

  virtual LogSphericalSolver* Clone() const {
    return new LogSphericalSolver(*this);
  }

  virtual void operator()(Field<R, 3, DGrid>& E, Field<R, 3, DGrid>& B,
                          Field<RJ, 3, DGrid>& J,
                          const apt::Grid<R, DGrid>& grid,
                          const mpi::CartComm& cart, int timestep,
                          R dt) const override;

  static constexpr int min_guard(int iteration_in_inverting_operator) {
    return 2 + iteration_in_inverting_operator;
  }
};

}  // namespace field

#endif
