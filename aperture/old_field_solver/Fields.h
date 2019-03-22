#ifndef _FIELDS_H_
#define _FIELDS_H_

#include "MultiArray.h"
#include "Grid.h"

/// This is the base class for fields living on a grid. It maintains a
/// \ref Grid object and caches its linear size. It also implements a
/// function to check whether two grid extents are the same.
class FieldBase {
 public:
  /// Default constructor, initialize an empty grid and zero grid size.
  FieldBase() : _grid(), _gridSize(0) {}

  /// Main constructor, initializes the grid according to a given \ref
  /// Grid object.
  FieldBase(const Grid& grid) : _grid(grid), _gridSize(grid.size()) {}

  /// Destructor. No need to destruct anything since we didn't
  /// allocate dynamic memory. However it is useful to declare this as
  /// virtual since we need to derive from this class.
  virtual ~FieldBase() {}

  // Accessor methods for everything
  const Grid& grid() const { return _grid; }
  int gridSize() const { return _gridSize; }
  Extent extent() const { return _grid.extent(); }

 protected:
  Grid _grid;  ///< Grid that this field lives on
  int _gridSize;    //!< Cache the grid size for easier retrieval

  void check_grid_extent(const Extent& ext1, const Extent& ext2);
};  // ----- end of class FieldBase -----

/// Class for a scalar field with one component.
template <typename T>
class ScalarField : public FieldBase {
 public:
  typedef T data_type;
  typedef Grid grid_type;
  typedef MultiArray<T> array_type;
  typedef ScalarField<T> self_type;

  // Constructors and destructor
  ScalarField();
  ScalarField(const grid_type& grid);
  ScalarField(const self_type& field);
  ScalarField(self_type&& field);
  virtual ~ScalarField();

  // Core functions
  void initialize();
  void assign(data_type value);
  void copyFrom(const self_type& field);
  self_type& operator = ( const self_type& field);
  self_type& operator = ( self_type&& field);

  void resize(const Grid& grid);

  // Arithmetic operations
  self_type& multiplyBy(data_type value);
  self_type& multiplyBy(const ScalarField<T>& field);
  self_type& addBy(data_type value);
  self_type& addBy(const ScalarField<T>& field);
  self_type& subtractBy(data_type value);
  self_type& subtractBy(const ScalarField<T>& field);

  // Index operator
  data_type& operator()(int x, int y = 0, int z = 0) { return _array(x, y, z); }
  const data_type& operator()(int x, int y = 0, int z = 0) const {
    return _array(x, y, z);
  }

  // Accessor methods
  array_type& data() { return _array; }
  const array_type& data() const { return _array; }
  data_type* ptr() { return _array.data(); }
  const data_type* ptr() const { return _array.data(); }

 private:
  array_type _array;
};  // ----- end of class ScalarField -----

template <typename T>
class VectorField : public FieldBase {
 public:
  typedef T data_type;
  typedef Grid grid_type;
  typedef MultiArray<T> array_type;
  typedef VectorField<T> self_type;

  /// Constructors and Destructor
  VectorField();
  VectorField(const grid_type& grid);
  VectorField(const self_type& field);
  VectorField(self_type&& field);
  virtual ~VectorField();

  self_type& operator=(const self_type& field);
  self_type& operator=( self_type&& field);

  /// Core functions
  void initialize();
  void assign(data_type value, int n);
  void assign(data_type value);
  void copyFrom(const self_type& field);

  void resize(const Grid& grid);

  /// Arithmetic operations
  self_type& multiplyBy(data_type value);
  self_type& multiplyBy(const ScalarField<T>& field);
  self_type& addBy(data_type value, int n);
  self_type& addBy(const VectorField<T>& field);
  self_type& subtractBy(data_type value, int n);
  self_type& subtractBy(const VectorField<T>& field);

  /// Index operator
  data_type& operator()(int n, int x, int y = 0, int z = 0) {
    return _array[n](x, y, z);
  }
  const data_type& operator()(int n, int x, int y = 0, int z = 0) const {
    return _array[n](x, y, z);
  }

  /// Accessor methods
  array_type& data(int n) { return _array[n]; }
  const array_type& data(int n) const { return _array[n]; }
  data_type* ptr(int n) { return _array[n].data(); }
  const data_type* ptr(int n) const { return _array[n].data(); }

 private:
  array_type _array[VECTOR_DIM];
};  // ----- end of class VectorField -----

inline Index GetStagProperty(FieldType type, int comp) {
  if (0 == comp) {
      return FieldType::ETYPE == type ? Index(1,0,0) : Index(0,1,1);
  } else if (1 == comp) {
      return FieldType::ETYPE == type ? Index(0,1,0) : Index(1,0,1);
  } else {
      return FieldType::ETYPE == type ? Index(0,0,1) : Index(1,1,0);
  }
}

inline int DirIncrement (int dir, const Extent& ext) {
  switch (dir) {
    case 0:
      return 1;
    case 1:
      return ext.width();
    case 2:
      return ext.width() * ext.height();
    default:
      return 1;
  }
}

#endif  // ----- #ifndef _FIELDS_H_  -----
