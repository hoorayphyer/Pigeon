#ifndef _MULTIARRAY_H_
#define _MULTIARRAY_H_

#include "Types.h"
#include <type_traits>
#include <algorithm>

///
/// @file   MultiArray.h
/// @author Alex Chen <fizban007@gmail.com>
/// @date   Thu Mar 19 14:42:25 2015
///
/// @brief This file defines the MultiArray class which is a
/// multi-dimensional array used in the code to store field values on
/// a grid.
///


/// The MultiArray class is a unified interface for 1D, 2D and 3D
/// arrays with proper index access and memory management.
template <typename T>
class MultiArray {
 public:
  typedef T data_type;
  typedef T* ptr_type;
  typedef MultiArray<T> self_type;

  /// Iterator for the MultiArray class, defines both the const and
  /// nonconst iterators at the same time.
  template <bool isConst>
  class const_nonconst_iterator {
   public:
    // Typedefs for easier writing and const/nonconst switch
    typedef MultiArray<T> array_type;
    typedef const_nonconst_iterator<isConst> self_type;
    typedef Index index_type;
    typedef T data_type;
    typedef typename std::conditional<isConst, const T*, T*>::type ptr_type;
    typedef typename std::conditional<isConst, const T&, T&>::type ref_type;
    typedef typename std::conditional<isConst, const array_type&,
                                      array_type&>::type arr_ref_type;
    typedef typename std::conditional<isConst, const array_type*,
                                      array_type*>::type arr_ptr_type;

    // Constructors and destructors
    const_nonconst_iterator(arr_ref_type array, index_type pos)
        : _array(array), _pos(pos) {}

    const_nonconst_iterator(arr_ref_type array, int x, int y, int z)
        : _array(array), _pos(x, y, z) {}

    const_nonconst_iterator(arr_ref_type array, int idx)
        : _array(array), _pos(idx, _array.extent()) {
      // std::cout << "Index is " << idx << std::endl;
    }

    /// Copy constructor.
    const_nonconst_iterator(const self_type& iter)
        : _array(iter._array), _pos(iter._pos) {}

    ~const_nonconst_iterator() {}

    /// Add-assign operator, takes in 3 integers
    self_type& operator+=(const Vec3<int>& extent) {
      _pos.x += extent.x;
      _pos.y += extent.y;
      _pos.z += extent.z;
      return (*this);
    }

    /// Minus-assign operator, takes in 3 integers
    self_type& operator-=(const Vec3<int>& extent) {
      _pos.x -= extent.x;
      _pos.y -= extent.y;
      _pos.z -= extent.z;
      return (*this);
    }

    /// Add operator, takes in 3 integers
    self_type operator+(const Vec3<int>& extent) const {
      self_type tmp(*this);
      tmp += extent;
      return tmp;
    }

    /// Minus operator, takes in 3 integers
    self_type operator-(const Vec3<int>& extent) const {
      self_type tmp(*this);
      tmp -= extent;
      return tmp;
    }

    /// Comparison operator, greater than. Iterators should be able to
    /// compare to each other regardless of constness.
    template <bool Const>
    bool operator>(const const_nonconst_iterator<Const>& iter) const {
      return (_pos.index(_array.extent()) >
              iter.pos().index(iter.array().extent()));
    }

    /// Comparison operator, less than. Iterators should be able to
    /// compare to each other regardless of constness.
    template <bool Const>
    bool operator<(const const_nonconst_iterator<Const>& iter) const {
      return (_pos.index(_array.extent()) <
              iter.pos().index(iter.array().extent()));
    }

    ////////////////////////////////////////////////////////////////////////////////
    //  Access operators
    ////////////////////////////////////////////////////////////////////////////////

    /// Pointer dereference, returns regular or const reference
    /// depending on whether the iterator is const
    ref_type operator*() const { return _array(_pos); }

    /// Index operator, returns regular or const reference
    /// depending on whether the iterator is const
    ref_type operator()(int x, int y = 0, int z = 0) const {
      return _array(_pos.x + x, _pos.y + y, _pos.z + z);
    }

    /// Index operator, returns regular or const reference
    /// depending on whether the iterator is const
    ref_type operator()(const Index& idx) const{
        return operator()(idx.x, idx.y, idx.z);
    }

    /// Linear index operator, returns regular or const reference
    /// depending on whether the iterator is const
    ref_type operator[](int n) const {
      return _array[_pos.index(_array.extent()) + n];
    }

    /// Returns a const reference to the array of this iterator.
    const array_type& array() const { return _array; }

    /// Returns a const reference to the position that this iterator points to.
    const index_type& pos() const { return _pos; }

   private:
    arr_ref_type _array;                ///< A reference to the underlying array
    index_type _pos;                    ///< Position that this iterator points to
  };  // ----- end of class const_nonconst_iterator

  typedef const_nonconst_iterator<false> iterator;
  typedef const_nonconst_iterator<true> const_iterator;

  ////////////////////////////////////////////////////////////////////////////////
  //  Constructors
  ////////////////////////////////////////////////////////////////////////////////

  /// Default constructor, initializes `_size` to zero and `_data` to
  /// `nullptr`.
  MultiArray() : _data(nullptr), _size(0) {
    _extent.width() = 0;
    _extent.height() = 1;
    _extent.depth() = 1;
  }

  /// Main constructor, initializes with given width, height, and
  /// depth of the array. Allocate memory in the initialization.
  explicit MultiArray(int width, int height = 1, int depth = 1)
      : _extent{width, height, depth} {
    _size = _extent.size();
    find_dim();

    _data = new T[_size];
  }

  /// Alternative main constructor, takes in an @ref Extent object and
  /// initializes an array of the corresponding extent.
  explicit MultiArray(const Extent& extent)
      : MultiArray(extent.width(), extent.height(), extent.depth()) {}

  /// Standard copy constructor.
  MultiArray(const self_type& other)
      : _extent(other._extent), _size(other._size), _dim(other._dim) {
    _data = new T[_size];
    copyFrom(other);
  }

  /// Standard move constructor.
  MultiArray(self_type&& other)
      : _extent(other._extent), _size(other._size), _dim(other._dim) {
    _data = other._data;
    other._data = nullptr;
    other._size = 0;
    other._extent.width() = 0;
    other._extent.height() = 1;
    other._extent.depth() = 1;
  }

  /// Destructor. Delete the member data array.
  virtual ~MultiArray() {
    if (_data != nullptr) {
      delete[] _data;
      _data = nullptr;
    }
  }

  /// Assignment operators for copying
  self_type& operator=(const self_type& other) {
    if (_extent != other._extent) {
      resize(other._extent);
    }
    copyFrom(other);
    return (*this);
  }

  /// Move assignment operator
  self_type& operator=(self_type&& other) {
    if (_extent != other._extent) {
      _extent = other._extent;
      _size = _extent.size();
      find_dim();
    }
    // If the memory is already allocated, then pointing _data to
    // another place will lead to memory leak.
    if (_data != nullptr) {
      delete[] _data;
    }
    _data = other._data;
    other._data = nullptr;
    return (*this);
  }

  /// Linearized indexing operator, read only
  const data_type& operator[](int idx) const {
    return _data[idx];
  }

  /// Linearized indexing operator, read and write
  data_type& operator[](int idx) {
    return _data[idx];
  }

  /// Vector indexing operator, read only
  const data_type& operator()(int x, int y = 0, int z = 0) const {
    int idx = x + y * _extent.width() + z * _extent.width() * _extent.height();
    return _data[idx];
  }

  /// Vector indexing operator, read and write
  data_type& operator()(int x, int y = 0, int z = 0) {
    int idx = x + y * _extent.width() + z * _extent.width() * _extent.height();
    return _data[idx];
  }

  /// Vector indexing operator using an @ref Index object, read only
  const data_type& operator()(const Index& index) const {
    int idx = index.index(_extent);
    return _data[idx];
  }

  /// Vector indexing operator using an @ref Index object, read and write
  data_type& operator()(const Index& index) {
    int idx = index.index(_extent);
    return _data[idx];
  }

  /// Copying the entire content from another vector
  void copyFrom(const self_type& other) {
    std::copy_n(other._data, _size, _data);
  }

  /// Set the whole array to a single initial value
  void assign(const data_type& value) {
    std::fill_n(_data, _size, value);
  }

  /// Resize the array.
  void resize(int width, int height = 1, int depth = 1) {
    _extent.width() = width;
    _extent.height() = height;
    _extent.depth() = depth;
    _size = _extent.size();
    find_dim();
    if (_data != nullptr) {
      delete[] _data;
    }
    _data = new T[_size];
    assign( static_cast<T>(0) );
  }

  /// Resize the array according to an \ref Extent object.
  void resize(Extent extent) {
    resize(extent.width(), extent.height(), extent.depth());
  }

  ////////////////////////////////////////////////////////////////////////////////
  //  Iterators
  ////////////////////////////////////////////////////////////////////////////////

  iterator begin() { return iterator(*this, 0); }
  iterator end() { return iterator(*this, _size); }
  iterator index(int x, int y = 0, int z = 0) {
    return iterator(*this, x, y, z);
  }
  iterator index(const Index& index) {
    return iterator(*this, index.x, index.y, index.z);
  }

  const_iterator begin() const { return const_iterator(*this, 0); }
  const_iterator end() const { return const_iterator(*this, _size); }
  const_iterator index(int x, int y = 0, int z = 0) const {
    return const_iterator(*this, x, y, z);
  }
  const_iterator index(const Index& index) const {
    return const_iterator(*this, index.x, index.y, index.z);
  }

  /// Get the dimensions of this array
  ///
  /// @return Dimension of the multi-array
  ///
  int dim() const { return _dim; }

  // Returns various sizes of the array
  int width() const { return _extent.width(); }
  int height() const { return _extent.height(); }
  int depth() const { return _extent.depth(); }
  int size() const { return _size; }
  const Extent& extent() const { return _extent; }

  /// Direct access to the encapsulated pointer
  T* data() { return _data; }
  const T* data() const { return _data; }

 private:
  void find_dim() {
    if (_extent.height() <= 1 && _extent.depth() <= 1)
      _dim = 1;
    else if (_extent.depth() <= 1)
      _dim = 2;
    else
      _dim = 3;
  }

  ptr_type _data;                       ///< Pointer to the data stored

  Extent _extent;                       ///< Extent of the array in all dimensions
  int _size;                            ///< Total size of the array
  int _dim;                             ///< Dimension of the array

};  // ----- end of class MultiArray -----

#endif  // ----- #ifndef _MULTIARRAY_H_  -----
