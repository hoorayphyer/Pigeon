#include "simulator/conf_file.hpp"

#include <cstdint>  // for int64
#include <format>

#include "toml++/toml.h"

using namespace std::string_view_literals;

namespace pic {

namespace {
using Node_t = toml::table;

template <typename T>
constexpr auto get_toml_native_type() {
  if constexpr (std::is_floating_point_v<T>) {
    return double{};
  } else if constexpr (std::is_integral_v<T>) {
    return int64{};
  } else {
    return T{};
  }
}

template <typename T>
using toml_native_t = std::remove_cvref_t<decltype(get_toml_native_type<T>())>;
}  // namespace

ConfFile ConfFile::load(const std::string& file) {
  ConfFile res;
  try {
    res.m_node = toml::parse_file(file);
  } catch (const toml::parse_error& err) {
    auto msg =
        std::format("Error parsing file '{}':\n{}\n  ({})"sv,
                    *err.source().path, err.description(), err.source().begin);

    throw std::runtime_error(msg);
  }

  return res;
}

ConfFile ConfFile::operator[](const std::string& entry) {
  ConfFile res;
  // c++20 has std::format for compile-time format and std::vformat for
  // runtime format
  res.m_current_entries = std::format("{}[{}]", m_current_entries, entry);
  res.m_node = std::any_cast<Node_t&>(m_node)[entry];
  return res;
}

template <typename T>
T ConfFile::as<T>() const {
  auto toml_native_val =
      std::any_cast<const Node_t&>(m_node).template value<toml_native_t<T>>();
  if (!toml_native_val) {
    auto msg =
        std::format("ERROR : setting {} with non-existent value in config file",
                    m_current_entries);
    throw std::runtime_error(msg);
  }

  return static_cast<T>(*toml_native_val);
}

template <typename T>
T ConfFile::as_or<T>(T val_default) const {
  auto toml_native_val = std::any_cast<const Node_t&>(m_node).template value_or(
      static_cast<toml_native_t<T>>(val_default));

  return static_cast<T>(toml_native_val);
}

#define INSTANTIATE_AS_TYPE(_TYPE_)                      \
  template _TYPE_ ConfFile::as_or<_TYPE_>(_TYPE_) const; \
  template _TYPE_ ConfFile::as<_TYPE_>() const

INSTANTIATE_AS_TYPE(bool);
INSTANTIATE_AS_TYPE(float);
INSTANTIATE_AS_TYPE(double);
INSTANTIATE_AS_TYPE(std::string);
INSTANTIATE_AS_TYPE(int);

}  // namespace pic
