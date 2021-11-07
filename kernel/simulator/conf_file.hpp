#pragma once
#include <any>
#include <string>

namespace pic {
struct ConfFile {
 public:
  static ConfFile load(const std::string& file);

  ConfFile operator[](const std::string& entry);

  template <typename T>
  T as() const;

  template <typename T>
  T as_or(T val_default) const;

 private:
  ConfFile() = default;
  std::string m_current_entries = "";
  std::any m_node;  // type erasure
  // std::any only stores copyable objects
};
}  // namespace pic
