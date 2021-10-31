#pragma once

#include <cstdio>
#include <string>

#include "apt/array.hpp"
#include "apt/pair.hpp"
#include "apt/vec_expression.hpp"

namespace std {
template <typename T, int N>
string to_string(const apt::array<T, N>& c) {
  string res = "( " + to_string(c[0]);
  for (int i = 1; i < N; ++i) res += (", " + to_string(c[i]));
  res += " )";
  return res;
}

template <typename T>
string to_string(const apt::pair<T>& c) {
  return "( " + to_string(c[0]) + ", " + to_string(c[1]) + " )";
}

template <typename E>
string to_string(const apt::VecExpression<E>& c) {
  string res = "( " + to_string(c[0]);
  for (int i = 1; i < E::NDim; ++i) res += (", " + to_string(c[i]));
  res += " )";
  return res;
}
}  // namespace std

template <typename OStream, typename T, int N>
OStream& operator<<(OStream& os, const apt::array<T, N>& c) {
  os << "( " << c[0];
  for (int i = 1; i < N; ++i) os << ", " << c[i];
  os << " )";
  return os;
}

template <typename OStream, typename T>
OStream& operator<<(OStream& os, const apt::pair<T>& c) {
  os << "( " << c[0] << ", " << c[1] << " )";
  return os;
}

template <typename OStream, typename E>
OStream& operator<<(OStream& os, const apt::VecExpression<E>& c) {
  os << "( " << c[0];
  for (int i = 1; i < E::NDim; ++i) os << ", " << c[i];
  os << " )";
  return os;
}

namespace apt {
template <typename... T>
std::string fmt(std::string s, T&&... x) {
  constexpr int SIZE = 20;
  char str[SIZE];
  std::snprintf(str, SIZE, s.c_str(), std::forward<T>(x)...);
  return {str};
}
}  // namespace apt
