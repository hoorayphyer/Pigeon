#pragma once
#include <any>
#include <optional>
#include <string>

namespace pic {
struct ConfFile {
 public:
  static ConfFile load(const std::string& file);

  ConfFile operator[](const std::string& entry) const;

  ConfFile operator[](int entry) const;

  template <typename T>
  T as() const;

  template <typename T>
  T as_or(T val_default) const;

  template <typename T>
  std::optional<T> optional() const;

 private:
  ConfFile() = default;
  std::string m_current_entries = "";
  std::any m_node;  // type erasure
  // std::any only stores copyable objects
};
}  // namespace pic
