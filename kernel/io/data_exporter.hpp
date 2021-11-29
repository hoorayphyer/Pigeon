#pragma once
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "apt/range.hpp"
#include "io/exportee.hpp"
#include "io/exportee_by_function.hpp"

namespace dye {
template <int>
struct Ensemble;
}

namespace io {

template <typename RealDS, int DGrid, typename Real,
          template <typename> class S, typename RealJ>
class DataExporter {
 public:
  DataExporter() = default;
  DataExporter(const DataExporter &other) = delete;
  DataExporter(DataExporter &&other) = default;

  auto &set_is_collinear_mesh(bool collinear) {
    m_is_collinear_mesh = collinear;
    return *this;
  }

  auto &set_downsample_ratio(int ratio) {
    m_downsample_ratio = ratio;
    return *this;
  }

  auto &set_data_dir(std::string dir) {
    m_data_dir = std::move(dir);
    return *this;
  }

  auto &set_mesh_name(std::string name) {
    m_mesh_name = std::move(name);
    return *this;
  }

  auto &set_range(apt::array<apt::Range, DGrid> range) {
    m_range = std::move(range);
    return *this;
  }

  auto &add_exportee(FieldExportee<RealDS, DGrid, Real, RealJ> *exportee) {
    m_fexps.emplace_back(exportee);
    return *this;
  }

  auto &add_exportee(PtcExportee<RealDS, DGrid, Real, S> *exportee) {
    m_pexps.emplace_back(exportee);
    return *this;
  }

  using FieldExporteeByFunc = FexpTbyFunction<RealDS, DGrid, Real, RealJ>;
  using ParticleExporteeByFunc = PexpTbyFunction<RealDS, DGrid, Real, S>;

  /*
   * overloads for field exportee by function
   */
  auto &add_exportee(std::string name, int num_comps,
                     typename FieldExporteeByFunc::F action,
                     typename FieldExporteeByFunc::H post_hook = nullptr) {
    m_fexps.emplace_back(
        new FieldExporteeByFunc(std::move(name), num_comps, action, post_hook));
    return *this;
  }

  /*
   * overloads for particle exportee by function
   */
  auto &add_exportee(std::string name, int num_comps,
                     typename ParticleExporteeByFunc::F action,
                     typename ParticleExporteeByFunc::H post_hook = nullptr) {
    m_pexps.emplace_back(new ParticleExporteeByFunc(std::move(name), num_comps,
                                                    action, post_hook));
    return *this;
  }

  auto &apply(std::function<void(DataExporter &)> f) {
    f(*this);
    return *this;
  }

  const auto &get_range() const { return m_range; }
  auto &get_range() { return m_range; }

  void export_data(
      int timestep, Real dt, int num_files,
      const std::optional<mpi::CartComm> &cart_opt,
      const dye::Ensemble<DGrid> &ens,
      const apt::Grid<Real, DGrid> &grid,  // local grid
      const field::Field<Real, 3, DGrid> &E,
      const field::Field<Real, 3, DGrid> &B,
      const field::Field<RealJ, 3, DGrid> &J,  // J is Jmesh on a replica
      const particle::map<particle::array<Real, S>> &particles,
      const particle::map<particle::Properties> &properties) const;

 private:
  bool m_is_collinear_mesh = false;
  int m_downsample_ratio = 1;
  const int m_mesh_ghost = 1;
  std::string m_data_dir;
  std::string m_mesh_name = "PICMesh";
  apt::array<apt::Range, DGrid> m_range;

  std::vector<std::unique_ptr<FieldExportee<RealDS, DGrid, Real, RealJ>>>
      m_fexps;
  std::vector<std::unique_ptr<PtcExportee<RealDS, DGrid, Real, S>>> m_pexps;
};
}  // namespace io
