#include <stdexcept>
#include "Fields.h"
#include "ArrayOperations.h"

using namespace Helpers;

void FieldBase::check_grid_extent(const Extent& ext1, const Extent& ext2) {
  if (ext1 != ext2) throw std::invalid_argument("Field grids don't match!");
}

////////////////////////////////////////////////////////////////////////////////
//  Scalar Field Implementation
////////////////////////////////////////////////////////////////////////////////

template <typename T>
ScalarField<T>::ScalarField()
    : FieldBase(), _array() {}

template <typename T>
ScalarField<T>::ScalarField(const grid_type& grid)
    : FieldBase(grid), _array(grid.extent()) {
  //  std::cout << grid.extent().depth() << std::endl;
}

template <typename T>
ScalarField<T>::ScalarField(const self_type& field)
    : FieldBase(field._grid), _array(field._array) {}

template <typename T>
ScalarField<T>::ScalarField(self_type&& field)
    : FieldBase(field._grid), _array(std::move(field._array)) {}

template <typename T>
ScalarField<T>::~ScalarField() {}

template <typename T>
void ScalarField<T>::initialize() {
  // Assign the field to zero, whatever 0 corresponds to for type #T
  _array.assign(static_cast<T>(0));
}

template <typename T>
void ScalarField<T>::assign(data_type value) {
  // Assign a uniform value to the array
  _array.assign(value);
}

template <typename T>
void ScalarField<T>::copyFrom(const self_type& field) {
  /// We can copy as long as the extents are the same
  this->check_grid_extent(this->_grid.extent(), field.grid().extent());

  _array.copyFrom(field._array);
}

template <typename T>
ScalarField<T>& ScalarField<T>::operator = ( const self_type& field){
    this->_grid = field._grid;
    this->_gridSize = field._gridSize;
    _array = field._array;
    return (*this);
}

template <typename T>
ScalarField<T>& ScalarField<T>::operator = ( self_type&& field){
    this->_grid = field._grid;
    this->_gridSize = field._gridSize;
    _array = std::move(field._array);
    return (*this);
}

template <typename T>
void ScalarField<T>::resize (const Grid& grid) {
  this->_grid = grid;
  this->_gridSize = grid.size();
  _array.resize(grid.extent());
}

template <typename T>
ScalarField<T>& ScalarField<T>::multiplyBy(data_type value) {
  map_multi_array(_array.begin(), this->_grid.extent(), Op_MultConst<T>(value));
  return (*this);
}

template <typename T>
ScalarField<T>& ScalarField<T>::multiplyBy(
    const ScalarField<T>& field) {
  this->check_grid_extent(this->_grid.extent(), field.grid().extent());

  map_multi_array(_array.begin(), field.data().begin(), this->_grid.extent(),
                  Op_MultAssign<T>());
  return (*this);
}

template <typename T>
ScalarField<T>& ScalarField<T>::addBy(data_type value) {
  map_multi_array(_array.begin(), this->_grid.extent(), Op_PlusConst<T>(value));
  return (*this);
}

template <typename T>
ScalarField<T>& ScalarField<T>::addBy(
    const ScalarField<T>& field) {
  this->check_grid_extent(this->_grid.extent(), field.grid().extent());

  map_multi_array(_array.begin(), field.data().begin(), this->_grid.extent(),
                  Op_PlusAssign<T>());
  return (*this);
}

template <typename T>
ScalarField<T>& ScalarField<T>::subtractBy(data_type value) {
  map_multi_array(_array.begin(), this->_grid.extent(), Op_MinusConst<T>(value));
  return (*this);
}

template <typename T>
ScalarField<T>& ScalarField<T>::subtractBy(const ScalarField<T> &field) {
  this->check_grid_extent(this->_grid.extent(), field.grid().extent());

  map_multi_array(_array.begin(), field.data().begin(), this->_grid.extent(),
                  Op_MinusAssign<T>());
  return (*this);

}

////////////////////////////////////////////////////////////////////////////////
//  Vector Field Implementation
////////////////////////////////////////////////////////////////////////////////

template <typename T>
VectorField<T>::VectorField()
    : FieldBase(), _array() {}

template <typename T>
VectorField<T>::VectorField(const grid_type& grid)
    : FieldBase(grid) {
  for (int i = 0; i < VECTOR_DIM; ++i) {
    _array[i] = std::move(array_type(grid.extent()));
  }
}

template <typename T>
VectorField<T>::VectorField(const self_type& field)
    : FieldBase(field._grid), _array(field._array) {}

template <typename T>
VectorField<T>::VectorField(self_type&& field)
    : FieldBase(field._grid), _array(std::move(field._array)) {}

template <typename T>
VectorField<T>::~VectorField() {}

template <typename T>
VectorField<T>& VectorField<T>::operator= (const self_type& other) {
    this->_grid = other._grid;
    this->_gridSize = other._gridSize;
    for (int i = 0; i < VECTOR_DIM; ++i)
      _array[i] = other._array[i];
    return (*this);
}

template <typename T>
VectorField<T>& VectorField<T>::operator= ( self_type&& other) {
    this->_grid = other._grid;
    this->_gridSize = other._gridSize;
    for (int i = 0; i < VECTOR_DIM; ++i)
      _array[i] = std::move(other._array[i]);
    return (*this);
}

template <typename T>
void VectorField<T>::initialize() {
  for (int i = 0; i < VECTOR_DIM; ++i) {
    _array[i].assign(static_cast<T>(0));
  }
}

template <typename T>
void VectorField<T>::assign(data_type value, int n) {
  _array[n].assign(value);
}

template <typename T>
void VectorField<T>::assign(data_type value) {
  for (int i = 0; i < VECTOR_DIM; i++) {
    _array[i].assign(value);
  }
}

template <typename T>
void VectorField<T>::copyFrom(const self_type& field) {
  /// We can copy as long as the extents are the same
  this->check_grid_extent(this->_grid.extent(), field.grid().extent());

  for (int i = 0; i < VECTOR_DIM; ++i) {
    _array[i].copyFrom(field._array[i]);
  }
}

template <typename T>
void VectorField<T>::resize (const Grid& grid) {
  this->_grid = grid;
  this->_gridSize = grid.size();
  for (int i = 0; i < VECTOR_DIM; i++) {
    _array[i].resize(grid.extent());
  }
}

template <typename T>
VectorField<T>& VectorField<T>::multiplyBy(data_type value) {
  for (int i = 0; i < VECTOR_DIM; ++i) {
    map_multi_array(_array[i].begin(), this->_grid.extent(),
                    Op_MultConst<T>(value));
  }
  return (*this);
}

template <typename T>
VectorField<T>& VectorField<T>::multiplyBy(
    const ScalarField<T>& field) {
  this->check_grid_extent(this->_grid.extent(), field.grid().extent());

  for (int i = 0; i < VECTOR_DIM; ++i) {
    map_multi_array(_array[i].begin(), field.data().begin(),
                    this->_grid.extent(), Op_MultAssign<T>());
  }
  return (*this);
}

template <typename T>
VectorField<T>& VectorField<T>::addBy(data_type value, int n) {
  map_multi_array(_array[n].begin(), this->_grid.extent(),
                  Op_PlusConst<T>(value));
  return (*this);
}

template <typename T>
VectorField<T>& VectorField<T>::addBy(
    const VectorField<T>& field) {
  this->check_grid_extent(this->_grid.extent(), field.grid().extent());

  for (int i = 0; i < VECTOR_DIM; ++i) {
    map_multi_array(_array[i].begin(), field.data(i).begin(),
                    this->_grid.extent(), Op_PlusAssign<T>());
  }
  return (*this);
}

template <typename T>
VectorField<T>& VectorField<T>::subtractBy(data_type value, int n) {
  map_multi_array(_array[n].begin(), this->_grid.extent(),
                  Op_MinusConst<T>(value));
  return (*this);
}

template <typename T>
VectorField<T>& VectorField<T>::subtractBy(const VectorField<T> &field) {
  this->check_grid_extent(this->_grid.extent(), field.grid().extent());

  for (int i = 0; i < VECTOR_DIM; ++i) {
    map_multi_array(_array[i].begin(), field.data(i).begin(),
                    this->_grid.extent(), Op_MinusAssign<T>());
  }
  return (*this);
}

////////////////////////////////////////////////////////////////////////////////
//  Explicit instantiations
////////////////////////////////////////////////////////////////////////////////

template class ScalarField<double>;
template class ScalarField<float>;
template class ScalarField<unsigned int>;

template class VectorField<double>;
template class VectorField<float>;
