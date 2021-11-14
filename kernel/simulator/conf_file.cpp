#include "simulator/conf_file.hpp"

#include <cassert>
#include <cstdint>  // for int64_t
#include <memory>
#include <sstream>  // for toml parsing error

#include "fmt/format.h"
#include "toml.hpp"

using namespace std::string_view_literals;

namespace pic {

namespace {

class Node_t {
 public:
  Node_t(toml::table&& table) { m_table.reset(new auto {std::move(table)}); }
  Node_t(const Node_t&) = default;
  Node_t(Node_t&&) = default;

  Node_t operator[](const std::string& entry) const {
    Node_t res;
    res.m_table = m_table;
    if (!m_node) {
      res.m_node.emplace((*m_table)[entry]);
    } else {
      res.m_node.emplace((*m_node)[entry]);
    }
    return res;
  }

  Node_t operator[](int entry) const {
    Node_t res;
    res.m_table = m_table;
    // m_table doesn't support accessing by an integer
    assert(m_node);
    res.m_node.emplace((*m_node)[entry]);
    return res;
  }

  auto& native() {
    assert(m_node);
    return *m_node;
  }

  const auto& native() const {
    assert(m_node);
    return *m_node;
  }

 private:
  Node_t() = default;
  using node_view_t = toml::node_view<toml::node>;
  std::shared_ptr<toml::table> m_table;
  std::optional<node_view_t> m_node;
};

template <typename T>
constexpr auto get_toml_native_type() {
  if constexpr (std::is_floating_point_v<T>) {
    return double{};
  } else if constexpr (std::is_integral_v<T>) {
    return std::int64_t{};
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
    res.m_node = Node_t{toml::parse_file(file)};
  } catch (const toml::parse_error& err) {
    // error.source().begin is type `source_positon`, which only has overload <<
    // for streams but not direct conversion to string. So use stringstream
    // here.
    std::ostringstream oss;
    oss << "Error parsing file '" << *err.source().path << "':\n"
        << err.description() << "\n  (" << err.source().begin << ")\n";

    throw std::runtime_error(oss.str());
  }

  return res;
}

ConfFile ConfFile::operator[](const std::string& entry) {
  ConfFile res;
  res.m_current_entries = fmt::format("{}[{}]", m_current_entries, entry);
  res.m_node = std::any_cast<Node_t&>(m_node)[entry];
  return res;
}

ConfFile ConfFile::operator[](int entry) {
  ConfFile res;
  res.m_current_entries = fmt::format("{}[{}]", m_current_entries, entry);
  res.m_node = std::any_cast<Node_t&>(m_node)[entry];
  return res;
}

template <typename T>
T ConfFile::as() const {
  auto toml_native_val = std::any_cast<const Node_t&>(m_node)
                             .native()
                             .template value<toml_native_t<T>>();
  if (!toml_native_val) {
    auto msg =
        fmt::format("ERROR : setting {} with non-existent value in config file",
                    m_current_entries);
    throw std::runtime_error(msg);
  }

  return static_cast<T>(*toml_native_val);
}

template <typename T>
T ConfFile::as_or(T val_default) const {
  auto toml_native_val =
      std::any_cast<const Node_t&>(m_node).native().template value_or(
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
