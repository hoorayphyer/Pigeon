#pragma once

namespace mpi {
struct CartComm;
}

namespace apt {
struct Range;
template <typename T, int D>
struct array;
}  // namespace apt

#include "field/field.hpp"

namespace field {
// Require sub_range <= field.range() <= domain_range. This allows the following
// 1. field's range == sub_range == domain_range
// 2. field's range == domain_range, but sync is called within an action range
// that is a subset
// 3. field's range == sub_range, which is smaller than the domain_range
template <typename T, int DField, int DGrid>
void copy_sync_guard_cells(Field<T, DField, DGrid>& field,
                           const mpi::CartComm& comm,
                           const apt::array<apt::Range, DGrid>& domain_range,
                           const apt::array<apt::Range, DGrid>& sub_range);

template <typename T, int DField, int DGrid>
inline void copy_sync_guard_cells(Field<T, DField, DGrid>& field,
                                  const mpi::CartComm& comm) {
  copy_sync_guard_cells(field, comm, field.mesh().range(),
                        field.mesh().range());
}

template <typename T, int DField, int DGrid>
void merge_sync_guard_cells(Field<T, DField, DGrid>& field,
                            const mpi::CartComm& comm,
                            const apt::array<apt::Range, DGrid>& domain_range,
                            const apt::array<apt::Range, DGrid>& sub_range);

template <typename T, int DField, int DGrid>
inline void merge_sync_guard_cells(Field<T, DField, DGrid>& field,
                                   const mpi::CartComm& comm) {
  merge_sync_guard_cells(field, comm, field.mesh().range(),
                         field.mesh().range());
}
}  // namespace field
