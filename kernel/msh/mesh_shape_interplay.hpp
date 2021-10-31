#pragma once

#include "apt/type_traits.hpp"
#include "apt/vec.hpp"  // TODOL solely for interpolate as return value
#include "field/field.hpp"

namespace msh {
// NOTE The standard grid of the bulk of the mesh is a rescaled and shifted
// version of the actual grid such that the spacing = 1 and the first cell has
// index = 0. This way, there is no need of grid information. "Standard" is
// borrowed from "standard normal distribution".
template <typename Grid, typename Vec>
inline auto to_standard(const Grid& grid, const Vec& qabs) noexcept {
  apt::array<apt::remove_cvref_t<typename Vec::element_type>, Vec::NDim>
      q_std;                   // NOTE DVec is used here
  apt::foreach<0, Grid::NDim>  // NOTE DGrid instead of DVec
      ([](auto& q, auto q_abs,
          const auto& g) noexcept { q = (q_abs - g.lower()) / g.delta(); },
       q_std, qabs, grid);
  apt::foreach<Grid::NDim, Vec::NDim>(
      [](auto& q, auto q_abs) noexcept { q = q_abs; }, q_std, qabs);
  return q_std;
}
}  // namespace msh

namespace msh {
// NOTE q_std refers to the same "standard" as above
template <typename T, int DGrid, int Dq, typename ShapeF>
T interpolate(const field::Component<T, DGrid>& fcomp,
              const apt::array<T, Dq>& q_std, const ShapeF& shapef) noexcept;

// NOTE q_std refers to the same "standard" as above
template <typename T, int DField, int DGrid, int Dq, typename ShapeF>
apt::Vec<T, DField> interpolate(const field::Field<T, DField, DGrid>& field,
                                const apt::array<T, Dq>& q_std,
                                const ShapeF& shapef) noexcept;

// the opposite of interpolate. `var` is +=ed to field // NOTE, this is not same
// as depositing current
template <typename T, int DField, int DGrid, int Dq, typename ShapeF>
void deposit(field::Field<T, DField, DGrid>& field, T frac,
             apt::array<T, DField> var, const apt::array<T, Dq>& q_std,
             const ShapeF& shapef) noexcept;
}  // namespace msh

// current deposition
namespace msh {
template <typename RealJ, int DField, int DGrid, typename ShapeF, typename U>
void deposit(field::Field<RealJ, DField, DGrid>& J, U charge_over_dt,
             const ShapeF& shapef,
             const apt::array<U, DField>& q0_std,  // NOTE its DField, not DGrid
             const apt::array<U, DField>& q1_std);

}
