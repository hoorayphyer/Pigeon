#pragma once
#include <string>
#include <vector>

#include "silopp/silo_optlist.hpp"

typedef struct DBfile DBfile;

namespace silo {
enum class MeshType : char { Rect = 0, Curv };
}

namespace silo {
template <typename file_t>
struct Operations {
 private:
  inline DBfile* _dbfile() noexcept {
    return static_cast<file_t&>(*this).operator DBfile*();
  }

 public:
  template <typename T>
  void put_mesh(std::string meshname, const std::vector<std::vector<T>>& coords,
                MeshType mt, OptList optlist = {});

  // FIXME ad hoc for LogSpherical2D
  template <typename T>
  void put_mesh_noncollinear(std::string meshname, const T* const coords[],
                             int* quadmesh_dims, int ndims,
                             OptList optlist = {});

  // put scalar field
  template <typename T>
  void put_var(std::string varname, std::string meshname, const T* vardata,
               const std::vector<int>& dims, OptList optlist = {});

  // put vector or tensor field
  template <typename T>
  void put_var(std::string varname, std::string meshname,
               const std::vector<const T*>& vardata,
               const std::vector<int>& dims, OptList optlist = {});

  void put_multimesh(std::string multimeshname, int nblock, std::string file_ns,
                     std::string block_ns, MeshType mt, OptList optlist = {});

  void put_multivar(std::string multivarname, int nblock, std::string file_ns,
                    std::string block_ns, OptList optlist = {});

  template <typename T>
  void write(std::string varname, const T* vardata,
             const std::vector<int>& dims);

  template <typename T>
  void write(std::string varname, const std::vector<T>& vardata);

  template <typename T>
  void write(std::string varname, T vardata);
};
}  // namespace silo
