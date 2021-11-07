#pragma once

#include <string>
#include <vector>

typedef struct DBfile DBfile;

namespace silo {
struct Slice {
  int begin{};
  int end{};
  int stride = 1;
};

template <typename file_t>
struct Reader {
 private:
  inline DBfile* _dbfile() noexcept {
    return static_cast<file_t&>(*this).operator DBfile*();
  }

 public:
  bool var_exists(std::string varname);
  int var_datatype(std::string varname);
  std::vector<int> var_dims(std::string varname);
  int var_length(
      std::string varname);  // TODO silo GetVarLength only returns int.

  void read(std::string varname, void* var);

  template <typename T>
  inline T read1(std::string varname) {
    T res = static_cast<T>(0);
    read(std::move(varname), &res);
    return res;
  }

  template <typename T>
  inline std::vector<T> read1d(std::string varname) {
    std::vector<T> res(var_length(varname));
    read(std::move(varname), res.data());
    res.shrink_to_fit();
    return res;
  }

  void readslice(std::string varname, const std::vector<Slice>& slice,
                 void* var);
};
}  // namespace silo
