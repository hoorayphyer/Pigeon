#ifndef _MULTIARRAYHELPERS_H_
#define _MULTIARRAYHELPERS_H_

#include "../MultiArray.h"

namespace Helpers {

////////////////////////////////////////////////////////////////////////////////
///  Mapping an operation over a multiarray.
////////////////////////////////////////////////////////////////////////////////
template <typename InputIterator, typename UnaryOp>
void map_multi_array(const InputIterator& it, const Extent& range, UnaryOp op) {
  for (int k = 0; k < range.depth(); ++k) {
    for (int j = 0; j < range.height(); ++j) {
      for (int i = 0; i < range.width(); ++i) {
        op(it(i, j, k));
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
///  Mapping an operation over two multiarrays.
////////////////////////////////////////////////////////////////////////////////
template <typename InputIterator, typename OutputIterator, typename BinaryOp>
void map_multi_array(const OutputIterator& output, const InputIterator& input,
                     const Extent& range, BinaryOp op) {
  for (int k = 0; k < range.depth(); ++k) {
    for (int j = 0; j < range.height(); ++j) {
      for (int i = 0; i < range.width(); ++i) {
        op(output(i, j, k), input(i, j, k));
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
/// class for element-wise multiplication of several arrays
/// (this class is more for fun)
////////////////////////////////////////////////////////////////////////////////
template <typename T1, typename T2>
class ArrayProduct {
 private:
  // These are arrays which implement parenthese operators with index
  // i, j, k or an Index object
  const T1* _array1;
  const T2* _array2;

 public:
  typedef ArrayProduct<T1, T2> self_type;

  template <bool isConst>
  class const_nonconst_iterator {
   public:
    typedef ArrayProduct<T1, T2> arrprd_type;
    typedef const_nonconst_iterator<isConst> self_type;
    typedef Index index_type;
    typedef Scalar data_type;

    typedef
        typename std::conditional<isConst, const data_type&, data_type&>::type
            data_ref_type;
    typedef
        typename std::conditional<isConst, const data_type*, data_type*>::type
            data_ptr_type;
    typedef typename std::conditional<isConst, const arrprd_type&,
                                      arrprd_type&>::type product_ref_type;
    typedef typename std::conditional<isConst, const arrprd_type*,
                                      arrprd_type*>::type product_ptr_type;

    /// Constructors and destructors
    const_nonconst_iterator(product_ref_type arrprd,
                            index_type pos)  // index ref for pos?????
        : _arrprd(arrprd), _pos(pos) {}
    const_nonconst_iterator(product_ref_type arrprd, int x, int y, int z)
        : _arrprd(arrprd), _pos(x, y, z) {}
    const_nonconst_iterator(product_ref_type arrprd, int index)
        : _arrprd(arrprd), _pos(index, _arrprd.extent()) {}
    const_nonconst_iterator(const self_type& iter)
        : _arrprd(iter._arrprd), _pos(iter._pos) {}
    ~const_nonconst_iterator() {}

    /// Addition operators
    self_type& operator+=(const Extent& extent) {
      _pos.x += extent.width();
      _pos.y += extent.height();
      _pos.z += extent.depth();
      return *this;
    }

    self_type& operator-=(const Extent& extent) {
      _pos.x -= extent.width();
      _pos.y -= extent.height();
      _pos.z -= extent.depth();
      return *this;
    }

    self_type operator+(const Extent& extent) const {
      self_type tmp(*this);
      tmp += extent;
      return tmp;
    }

    self_type operator-(const Extent& extent) const {
      self_type tmp(*this);
      tmp -= extent;
      return tmp;
    }

    /// Conditional operators
    template <bool Const>
    bool operator>(const_nonconst_iterator<Const>& iter) const {
      return (_pos.index(_arrprd._array1->extent()) >
              iter._pos.index(iter._arrprd._array1->extent()));
    }

    template <bool Const>
    bool operator<(const_nonconst_iterator<Const>& iter) const {
      return (_pos.index(_arrprd._array1->extent()) <
              iter._pos.index(iter._arrprd._array1->extent()));
    }

    /// Access operators
    data_type operator*() const { return _arrprd(_pos); }
    data_type operator()(int x, int y = 0, int z = 0) const {
      return _arrprd(_pos.x + x, _pos.y + y, _pos.z + z);
    }
    data_type operator()(const Index& idx) const {
      return operator()(idx.x, idx.y, idx.z);
    }

    const product_ref_type arrprd() const { return _arrprd; }
    const index_type& pos() const { return _pos; }

   private:
    product_ref_type _arrprd;
    index_type _pos;
  };  // ----- end of class const_nonconst_iterator

  typedef const_nonconst_iterator<false> iterator;
  typedef const_nonconst_iterator<true> const_iterator;

  /// Constructors and destructors
  ArrayProduct() : _array1(nullptr), _array2(nullptr) {}
  ArrayProduct(const T1& a1, const T2& a2) {
    if (sizeof(T1) > 2 * sizeof(T1*))
      _array1 = &a1;
    else
      _array1 = new T1(a1);

    if (sizeof(T2) > 2 * sizeof(T2*))
      _array2 = &a2;
    else
      _array2 = new T2(a2);
  }
  ArrayProduct(const self_type& other) {
    if (sizeof(T1) > 2 * sizeof(T1*))
      _array1 = other._array1;
    else
      _array1 = new T1(*(other._array1));

    if (sizeof(T2) > 2 * sizeof(T2*))
      _array2 = other._array2;
    else
      _array2 = new T2(*(other._array2));
  }

  ~ArrayProduct() {
    if (sizeof(T1) <= 2 * sizeof(T1*)) delete _array1;
    if (sizeof(T2) <= 2 * sizeof(T2*)) delete _array2;
  }

  /// Iterators
  iterator begin() { return iterator(*this, 0); }
  iterator end() { return iterator(*this, this->size()); }
  iterator index(int x, int y = 0, int z = 0) {
    return iterator(*this, x, y, z);
  }
  iterator index(const Index& index) {
    return iterator(*this, index.x, index.y, index.z);
  }
  const_iterator begin() const { return const_iterator(*this, 0); }
  const_iterator end() const { return const_iterator(*this, this->size()); }
  const_iterator index(int x, int y = 0, int z = 0) const {
    return const_iterator(*this, x, y, z);
  }
  const_iterator index(const Index& index) const {
    return const_iterator(*this, index.x, index.y, index.z);
  }

  /// Access operators
  inline Scalar operator()(int i, int j, int k) const {
    return (*_array1)(i, j, k) * (*_array2)(i, j, k);
  }
  inline Scalar operator()(Index idx) const {
    return (*_array1)(idx) * (*_array2)(idx);
  }
  inline Scalar operator[](int n) const {
    return (*_array1)[n] * (*_array2)[n];
  }

  /// Get the dimensions
  const Extent& extent() const { return _array1->extent(); }
  int size() const { return _array1->size(); }
  // ArrayGlue<T1, T2>& operator= (ArrayGlue<T1, T2>)
};  // ----- end of class ArrayGlue

////////////////////////////////////////////////////////////////////////////////
///  Generic operator * to glue two arrays together. Has several
///  overloaded cases.
////////////////////////////////////////////////////////////////////////////////
template <typename Array, typename T1, typename T2>
ArrayProduct<Array, ArrayProduct<T1, T2> > operator*(
    const Array& left, const ArrayProduct<T1, T2>& right) {
  return ArrayProduct<Array, ArrayProduct<T1, T2> >(left, right);
}

template <typename Array, typename T1, typename T2>
ArrayProduct<ArrayProduct<T1, T2>, Array> operator*(
    const ArrayProduct<T1, T2>& left, const Array& right) {
  return ArrayProduct<ArrayProduct<T1, T2>, Array>(left, right);
}

template <typename T>
ArrayProduct<MultiArray<T>, MultiArray<T> > operator*(
    const MultiArray<T>& left, const MultiArray<T>& right) {
  return ArrayProduct<MultiArray<T>, MultiArray<T> >(left, right);
}
}  // namespace Helpers

#endif  // ----- #ifndef _MULTIARRAYHELPERS_H_  -----
