#pragma once

#include <vector>

#include "apt/type_traits.hpp"
#include "field/mesh.hpp"
#include "field/offset.hpp"

namespace field {
template <typename T, int DGrid, bool Const = true>
struct Component {
 private:
  using vector_type = apt::cond_t<Const, const std::vector<T>, std::vector<T>>;
  vector_type& _data;
  const Mesh<DGrid>& _mesh;
  const apt::array<offset_t, DGrid>& _offset;

 public:
  constexpr Component(vector_type& data, const Mesh<DGrid>& mesh,
                      const apt::array<offset_t, DGrid>& offset) noexcept
      : _data(data), _mesh(mesh), _offset(offset) {}

  // allow conversion from non-const to const
  template <typename U = Component<T, DGrid, not Const>>
  constexpr Component(const std::enable_if_t<Const, U>& c) noexcept
      : Component(c.data(), c.mesh(), c.offset()) {}

  inline T& operator()(const apt::Index<DGrid>& i_bulk) {
    return _data[_mesh.linear_index(i_bulk)];
  }

  inline const T& operator()(const apt::Index<DGrid>& i_bulk) const {
    return _data[_mesh.linear_index(i_bulk)];
  }

  inline T& operator[](int i) { return _data[i]; }

  inline const T& operator[](int i) const { return _data[i]; }

  inline const auto& data() const noexcept { return _data; }
  inline auto& data() noexcept { return _data; }

  inline const auto& offset() const noexcept { return _offset; }
  inline const auto& mesh() const noexcept { return _mesh; }
};
}  // namespace field

namespace ckpt {
template <typename T, int DGrid>
struct FieldCkpt;
}

namespace field {

template <typename T, int DField, int DGrid>
struct Field {
 private:
  apt::array<std::vector<T>, DField> _comps;
  Mesh<DGrid> _mesh;
  apt::array<apt::array<offset_t, DGrid>, DField> _offset;

 public:
  using element_type = T;
  static constexpr int NDim = DField;
  friend class ckpt::FieldCkpt<T, DGrid>;

  Field() = default;
  // Field( const Field& ) = default;
  // Field( Field&& ) noexcept = default;

  Field(const Mesh<DGrid>& mesh) { resize(mesh); }

  inline Field& set_offset(int component,
                           const apt::array<offset_t, DGrid>& offset) noexcept {
    _offset[component] = offset;
    return *this;
  }

  inline Field& set_offset(int component, int ith_dim,
                           offset_t ofs_v) noexcept {
    _offset[component][ith_dim] = ofs_v;
    return *this;
  }

  inline const auto& mesh() const noexcept { return _mesh; }

  inline const auto operator[](int i) const noexcept {
    return Component<T, DGrid>(_comps[i], _mesh, _offset[i]);
  }

  inline auto operator[](int i) noexcept {
    return Component<T, DGrid, false>(_comps[i], _mesh, _offset[i]);
  }

  inline void reset() {
    for (int C = 0; C < DField; ++C)
      std::fill(_comps[C].begin(), _comps[C].end(), 0.0);
  }

  inline void resize(const Mesh<DGrid>& mesh) noexcept {
    auto size = mesh.stride().back();
    for (int C = 0; C < DField; ++C) {
      _comps[C].reserve(size);
      _comps[C].resize(size);
    }
    _mesh = mesh;
    reset();
  }
};
}  // namespace field
