#include <algorithm>

#include "filesys/filesys.hpp"
#include "silopp/silo++.hpp"

int main(int argc, char* argv[]) {
  if (argc != 3) {
    return 1;
  }
  std::string dir(argv[1]);
  std::string output(argv[2]);

  std::vector<int> supergrid = {4096, 4096};
  std::vector<int> guard = {2, 2};
  std::vector<bool> periodic = {false, false};
  std::vector<int> topos;
  std::vector<int> subgrid;
  std::vector<bool> ens_visited;
  int count_ens_unvisited = 0;

  std::array<std::vector<float>, 3> E, B;
  std::size_t size = 1;
  for (int i = 0; i < supergrid.size(); ++i)
    size *= supergrid[i] + 2 * guard[i];
  for (int i = 0; i < 3; ++i) {
    E[i].resize(size, {});
    B[i].resize(size, {});
    E[i].shrink_to_fit();
    B[i].shrink_to_fit();
  }

  for (auto f : fs::directory_iterator(dir)) {
    auto sf = silo::open(f, silo::Mode::Read);
    if (topos.size() == 0) {
      int ndims = sf.var_length("/cartesian_topology");
      topos.resize(ndims);
      sf.read("/cartesian_topology", topos.data());
      subgrid.resize(ndims);
      int N = 1;
      for (int i = 0; i < ndims; ++i) {
        subgrid[i] = supergrid[i] / topos[i];
        N *= topos[i];
      }
      count_ens_unvisited = N;
      ens_visited.resize(N, false);
    }

    auto load_subfield = [&guard, &subgrid, &supergrid](
                             auto& sf, const auto& coords, auto comp,
                             auto& superf) {
      std::vector<float> subf(sf.var_length(comp));
      sf.read(comp, subf.data());

      for (int j = -guard[1]; j < subgrid[1] + guard[1]; ++j) {
        int li = (j + guard[1]) * (subgrid[0] + 2 * guard[0]);
        int Li =
            (coords[0] * subgrid[0]) + (coords[1] * subgrid[1] + j + guard[1]) *
                                           (supergrid[0] + 2 * guard[0]);
        int ext = subgrid[0] + 2 * guard[0];
        std::copy_n(subf.data() + li, ext, superf.data() + Li);
      }
    };

    for (const auto& dname : sf.toc_dir()) {
      if (dname.find("ensemble") != 0) continue;

      int l = sf.read1<int>(dname + "/label");
      if (ens_visited[l]) continue;

      std::vector<int> coords(topos.size());
      sf.read(dname + "/cartesian_coordinates", coords.data());

      sf.cd(dname);
      if (sf.var_exists("rank0")) {
        for (int c = 0; c < 3; ++c) {
          load_subfield(sf, coords, "E" + std::to_string(c + 1), E[c]);
          load_subfield(sf, coords, "B" + std::to_string(c + 1), B[c]);
        }

        ens_visited[l] = true;
        if (--count_ens_unvisited == 0) break;
      }
      sf.cd("..");
    }
  }

  auto sf = silo::open(output, silo::Mode::Write);

  sf.write("supergrid", supergrid);
  sf.write("guard", guard);
  sf.write("E1", E[0]);
  sf.write("E2", E[1]);
  sf.write("E3", E[2]);

  sf.write("B1", B[0]);
  sf.write("B2", B[1]);
  sf.write("B3", B[2]);

  return 0;
}
