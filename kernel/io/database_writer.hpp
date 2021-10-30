#pragma once
#include <optional>
#include <string>

#include "filesys/filesys.hpp"
#include "mpipp/mpi++.hpp"
#include "silopp/pmpio.hpp"
#include "silopp/silo++.hpp"

namespace io {
class DatabaseWriter {
 private:
  const std::optional<mpi::CartComm>& m_cart_opt;

  // following are significant only on cart.rank() == 0
  std::string m_file_ns;
  std::string m_block_ns_prefix;
  std::string m_block_ns_postfix;

  silo::Pmpio m_pmpio;  // only significant on carts
  std::unique_ptr<silo::file_t>
      m_master{};  // only significant on cart.rank() == 0, use reference
                   // because silo library is sloppy about const-ness.
  std::string m_meshname;
  silo::OptList m_optlist;

  void set_namescheme(std::string str_ts, int num_files) {
    if (!m_cart_opt) return;

    // set up file_ns and block_ns. n below is thought of as the cartesian rank.
    // Only significant on world.rank() == 0
    constexpr char delimiter = '|';
    // NOTE file namescheme uses relative path so that when the directory is
    // moved, the data is still viewable
    {
      m_file_ns = delimiter + std::string("data/timestep") + str_ts +
                  "/set%d.silo" + delimiter + "n%" + std::to_string(num_files);
    }

    {
      m_block_ns_prefix = delimiter + std::string("cart");
      auto [c, topos] = m_cart_opt->coords_topos();
      for (int i = 0; i < topos.size(); ++i) m_block_ns_prefix += "_%03d";
      m_block_ns_prefix += "/";

      m_block_ns_postfix = "";
      // mpi uses row major numbering to map cartesian rank to coordinates, e.g.
      // (0,0) -> 0, (0,1) -> 1, (1,0) -> 2, (1,1) -> 3
      std::vector<int> strides(topos.size() +
                               1);  // strides = ( DzDyDx, DzDy, Dz, 1 )
      strides.back() = 1;
      for (int i = topos.size() - 1; i > -1; --i)
        strides[i] = strides[i + 1] * topos[i].dim();

      for (int i = 0; i < topos.size(); ++i) {
        m_block_ns_postfix +=
            delimiter + std::string("(n%" + std::to_string(strides[i]) + ")/") +
            std::to_string(strides[i + 1]);
      }
    }
  }

  template <typename T, int DGrid, bool C>
  void write_var_impl(std::string varname,
                      const field::Component<T, DGrid, C>& fcomp) const {
    if (!m_cart_opt) return;
    std::vector<int> dims(DGrid);
    for (int i = 0; i < DGrid; ++i) dims[i] = fcomp.mesh().range(i).full_size();
    m_pmpio([this, varname, &fcomp, dims = std::move(dims)](auto& dbfile) {
      dbfile.put_var(varname, m_meshname, fcomp.data().data(), dims);
    });

    if (m_master)
      m_master->put_multivar(varname, m_cart_opt->size(), m_file_ns,
                             m_block_ns_prefix + varname + m_block_ns_postfix,
                             m_optlist);
  }

 public:
  DatabaseWriter(const std::optional<mpi::CartComm>& cart_opt,
                 const std::string& prefix, int timestep, double dt,
                 int num_files)
      : m_cart_opt(cart_opt) {
    char str_ts[10];
    snprintf(str_ts, 10, "%07d", timestep);

    {
      if (m_cart_opt) {
        m_meshname = "PICMesh";

        m_pmpio.filename = prefix + "/data/timestep" + str_ts + "/set" +
                           std::to_string(m_cart_opt->rank() % num_files) +
                           ".silo";
        m_pmpio.dirname = "cart";
        for (const auto& x : m_cart_opt->coords()) {
          char tmp[10];
          sprintf(tmp, "_%03d", x);
          m_pmpio.dirname += tmp;
        }
        fs::mpido(*m_cart_opt, [&]() {
          fs::create_directories(prefix + "/data/timestep" + str_ts);
          m_master.reset(new silo::file_t(silo::open(
              prefix + "/timestep" + str_ts + ".silo", silo::Mode::Write)));
          set_namescheme(str_ts, num_files);
        });
        m_pmpio.comm = m_cart_opt->split((m_cart_opt->rank()) % num_files);
        if (m_pmpio.comm and m_pmpio.comm->rank() == 0) {
          auto dbfile = open(m_pmpio.filename, m_pmpio.mode);
          const auto [cc, topos] = m_cart_opt->coords_topos();
          std::vector<int> tp(topos.size());
          for (int i = 0; i < tp.size(); ++i) tp[i] = topos[i].signed_dim();
          dbfile.write("/cartesian_topology", tp);
        }

        m_optlist[silo::Opt::TIME] = timestep * dt;
        m_optlist[silo::Opt::CYCLE] = timestep;
      }
    }
  }

  template <typename RealDS, typename Real, int DGrid>
  void write_mesh(const field::Mesh<DGrid>& mesh, bool is_collinear_mesh,
                  int mesh_ghost, int downsample_ratio,
                  const apt::Grid<Real, DGrid>& grid) const {
    if (!m_cart_opt) return;

    // set quadmesh ghost cells
    /* this is the most correct way to do ghost, but it needs mesh to support
    variable guards
    // std::vector<int> lo_ofs( DGrid, 0);
    // std::vector<int> hi_ofs( DGrid, 0);
    // {
    //   const auto& c = ens.cart_coords;
    //   for ( int i = 0; i < DGrid; ++i ) {
    //     if ( c[i] > 0 ) lo_ofs[i] = silo_mesh_ghost;
    //     if ( c[i] < ens.cart_dims[i] - 1 ) hi_ofs[i] = silo_mesh_ghost;
    //   }
    // }
    */
    std::vector<int> lo_ofs(DGrid, mesh_ghost);
    std::vector<int> hi_ofs(DGrid, mesh_ghost);
    auto optlist_mesh = m_optlist;
    optlist_mesh[silo::Opt::LO_OFFSET] = lo_ofs;
    optlist_mesh[silo::Opt::HI_OFFSET] = hi_ofs;
    optlist_mesh[silo::Opt::BASEINDEX] =
        m_cart_opt->coords();  // need this in rectilinear mesh

    int quadmesh_dims[DGrid] = {};
    for (int i = 0; i < DGrid; ++i)
      quadmesh_dims[i] = mesh.range(i).size() / downsample_ratio + 1 +
                         lo_ofs[i] +
                         hi_ofs[i];  // plus 1 to include the upper boundary

    if (is_collinear_mesh) {
      std::vector<std::vector<RealDS> > coords(DGrid);
      for (int i = 0; i < DGrid; ++i) {
        auto& c = coords[i];
        auto dim = quadmesh_dims[i];
        c.reserve(dim);
        c.resize(dim, {});

        for (int j = 0; j < dim; ++j) c[j] = grid[i].absc(downsample_ratio * j);
      }

      m_pmpio([&](auto& dbfile) {
        dbfile.put_mesh(m_meshname, coords, silo::MeshType::Rect, optlist_mesh);
      });
    } else {
      // TODO these are solely for LogSpherical2D. It is a hotfix on visit
      // operators problems
      static_assert(DGrid == 2);

      RealDS* coords[DGrid];
      coords[0] = new RealDS[quadmesh_dims[0] * quadmesh_dims[1]];
      coords[1] = new RealDS[quadmesh_dims[0] * quadmesh_dims[1]];

      for (int j = 0; j < quadmesh_dims[1]; ++j) {
        for (int i = 0; i < quadmesh_dims[0]; ++i) {
          auto r = std::exp(grid[0].absc(downsample_ratio * (i - lo_ofs[0])));
          auto theta = grid[1].absc(downsample_ratio * (j - lo_ofs[1]));
          coords[0][i + j * quadmesh_dims[0]] = r * std::sin(theta);
          coords[1][i + j * quadmesh_dims[0]] = r * std::cos(theta);
        }
      }

      m_pmpio([&](auto& dbfile) {
        dbfile.put_mesh_noncollinear(m_meshname, coords, quadmesh_dims, DGrid,
                                     optlist_mesh);
      });

      delete[] coords[0];
      delete[] coords[1];
    }

    auto mesh_type =
        is_collinear_mesh ? silo::MeshType::Rect : silo::MeshType::Curv;
    if (m_master)
      m_master->put_multimesh(
          m_meshname, m_cart_opt->size(), m_file_ns,
          m_block_ns_prefix + m_meshname + m_block_ns_postfix, mesh_type,
          m_optlist);  // NOTE use m_optlist here, not optlist_mesh
  }

  template <typename RealDS, int DGrid>
  void write_var(const std::string& name, const std::string& name_postfix,
                 int field_dim,
                 const field::Field<RealDS, 3, DGrid>& fds) const {
    for (int comp = 0; comp < field_dim; ++comp) {
      auto varname = name;
      if (field_dim > 1) varname += std::to_string(comp + 1);
      if (!name_postfix.empty()) varname += "_" + name_postfix;
      write_var_impl(varname, fds[comp]);
    }
  }
};

}  // namespace io
