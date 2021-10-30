#ifndef _ARRAYOPERATIONS_H_
#define _ARRAYOPERATIONS_H_

// This file defines a series of operations on a chunk of a
// multi-array

//#include <stdexcept>
//#include "MultiArray.h"
#include <iostream>

#include "multiarray/MultiArrayHelpers.h"
#include "multiarray/Op.h"

template <typename It>
void check_bounds(const It& iterator, const Extent& extent) {
  // Check bounds
  if (iterator + extent - Extent(1, 1, 1) > iterator.array().end()) {
    std::cerr << iterator.pos() << std::endl;
    throw std::invalid_argument("Index out of bounds in array operation!");
  }
}

////////////////////////////////////////////////////////////////////////////////
///  Copy a chunk of data from input to output, with size specified
///  with extent.
////////////////////////////////////////////////////////////////////////////////
template <typename InputIt, typename OutputIt>
void copy(const OutputIt& output, const InputIt& input, const Extent& extent) {
  check_bounds(input, extent);
  check_bounds(output, extent);

  Helpers::map_multi_array(output, input, extent,
                           Op_Assign<typename InputIt::data_type>());
}

////////////////////////////////////////////////////////////////////////////////
///  Copy a chunk of data from input to output, with size specified
///  with extent, assuming the output is a linearized array.
////////////////////////////////////////////////////////////////////////////////
template <typename InputIt, typename OutputIt>
void copy_to_linear(const OutputIt& output, const InputIt& input,
                    const Extent& extent) {
  check_bounds(input, extent);
  // check_bounds(output, extent);

  // Helpers::map_multi_array(output, input, extent, Op_Assign<typename
  // InputIt::data_type>());
  for (int k = 0; k < extent.depth(); ++k) {
    for (int j = 0; j < extent.height(); ++j) {
      for (int i = 0; i < extent.width(); ++i) {
        int index = Index(i, j, k).index(extent);
        output[index] = input(i, j, k);
        // op(output(i, j, k), input(i, j, k));
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
///  Copy a chunk of data from input to output, with size specified
///  with extent, assuming the iutput is a linearized array.
////////////////////////////////////////////////////////////////////////////////
template <typename InputIt, typename OutputIt>
void copy_from_linear(const OutputIt& output, const InputIt& input,
                      const Extent& extent) {
  // check_bounds(input, extent);
  check_bounds(output, extent);

  // Helpers::map_multi_array(output, input, extent, Op_Assign<typename
  // InputIt::data_type>());
  for (int k = 0; k < extent.depth(); ++k) {
    for (int j = 0; j < extent.height(); ++j) {
      for (int i = 0; i < extent.width(); ++i) {
        int index = Index(i, j, k).index(extent);
        output(i, j, k) = input[index];
        // op(output(i, j, k), input(i, j, k));
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
///  Add a chunk of data from input to output, with size specified
///  with extent, assuming the iutput is a linearized array.
////////////////////////////////////////////////////////////////////////////////
template <typename InputIt, typename OutputIt>
void add_from_linear(const OutputIt& output, const InputIt& input,
                     const Extent& extent) {
  // check_bounds(input, extent);
  check_bounds(output, extent);

  // Helpers::map_multi_array(output, input, extent, Op_Assign<typename
  // InputIt::data_type>());
  for (int k = 0; k < extent.depth(); ++k) {
    for (int j = 0; j < extent.height(); ++j) {
      for (int i = 0; i < extent.width(); ++i) {
        int index = Index(i, j, k).index(extent);
        output(i, j, k) += input[index];
        // op(output(i, j, k), input(i, j, k));
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
///  Fill a chunk of input array with uniform value, the size of the
///  chunk given by extent.
////////////////////////////////////////////////////////////////////////////////
template <typename InputIt, typename T>
void fill(const InputIt& input, const Extent& extent, T value) {
  check_bounds(input, extent);

  Helpers::map_multi_array(input, extent,
                           Op_AssignConst<typename InputIt::data_type>(value));
}

////////////////////////////////////////////////////////////////////////////////
///  Multiply a chunk of the input array against the output array term
///  by term.
////////////////////////////////////////////////////////////////////////////////
template <typename InputIt, typename OutputIt>
void multiply(const OutputIt& output, const InputIt& input,
              const Extent& extent) {
  check_bounds(input, extent);
  check_bounds(output, extent);

  Helpers::map_multi_array(output, input, extent,
                           Op_MultAssign<typename InputIt::data_type>());
}

////////////////////////////////////////////////////////////////////////////////
///  Multiply a chunk of the input array by a single value.
////////////////////////////////////////////////////////////////////////////////
template <typename InputIt, typename T>
void multiply(const InputIt& input, const Extent& extent, T value) {
  check_bounds(input, extent);

  Helpers::map_multi_array(input, extent,
                           Op_MultConst<typename InputIt::data_type>(value));
}

////////////////////////////////////////////////////////////////////////////////
///  Add a chunk of the input array to the output array, size of the
///  chunk given by extent.
////////////////////////////////////////////////////////////////////////////////
template <typename InputIt, typename OutputIt>
void add(const OutputIt& output, const InputIt& input, const Extent& extent) {
  check_bounds(input, extent);
  check_bounds(output, extent);

  Helpers::map_multi_array(output, input, extent,
                           Op_PlusAssign<typename InputIt::data_type>());
}

////////////////////////////////////////////////////////////////////////////////
///  Add a chunk of the input array by a single value.
////////////////////////////////////////////////////////////////////////////////
template <typename InputIt, typename T>
void add(const InputIt& input, const Extent& extent, T value) {
  check_bounds(input, extent);

  Helpers::map_multi_array(input, extent,
                           Op_PlusConst<typename InputIt::data_type>(value));
}

////////////////////////////////////////////////////////////////////////////////
///  Subtract a chunk of the input array from the output array, size
///  of the chunk given by extent.
////////////////////////////////////////////////////////////////////////////////
template <typename InputIt, typename OutputIt>
void subtract(const OutputIt& output, const InputIt& input,
              const Extent& extent) {
  check_bounds(input, extent);
  check_bounds(output, extent);

  Helpers::map_multi_array(output, input, extent,
                           Op_MinusAssign<typename InputIt::data_type>());
}

////////////////////////////////////////////////////////////////////////////////
///  Subtract a chunk of the input array by a single value.
////////////////////////////////////////////////////////////////////////////////
template <typename InputIt, typename T>
void subtract(const InputIt& input, const Extent& extent, T value) {
  check_bounds(input, extent);

  Helpers::map_multi_array(input, extent,
                           Op_MinusConst<typename InputIt::data_type>(value));
}

#endif  // ----- #ifndef _ARRAYOPERATIONS_H_  -----
