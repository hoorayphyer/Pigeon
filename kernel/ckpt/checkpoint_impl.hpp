#include <cassert>

#include "ckpt/checkpoint.hpp"
#include "ckpt/checkpoint_helper.hpp"
#include "dye/dynamic_balance.hpp"
#include "field/field.hpp"
#include "filesys/filesys.hpp"
#include "mpipp/mpi++.hpp"
#include "particle/properties.hpp"
#include "silopp/pmpio.hpp"

#if PIC_DEBUG
#include "debug/debugger.hpp"
#include "logger/logger.hpp"
#endif

namespace ckpt {
// TODO save random seed
template <int DGrid, typename R, template <typename> class S>
std::string save_checkpoint(
    std::string prefix, const int num_parts,
    const std::optional<dye::Ensemble<DGrid>>& ens_opt, int timestep,
    const field::Field<R, 3, DGrid>& E, const field::Field<R, 3, DGrid>& B,
    const particle::map<particle::array<R, S>>& particles,
    const particle::map<particle::Properties>& properties) {
  bool is_idle = !ens_opt;
  int key = 0;
  if (!is_idle) key = ens_opt->label();
  // use ensemble label to put processes from same ensemble together
  auto active = mpi::world.split({is_idle}, key);
  if (is_idle) return "";
  const auto& ens = *ens_opt;

  char str_ts[10];
  sprintf(str_ts, "%06d\0", timestep);

  prefix = prefix + "/checkpoints/timestep" + str_ts;

  silo::Pmpio pmpio;

  int num_per_part =
      active->size() / num_parts + (active->size() % num_parts != 0);

  {
    pmpio.filename = prefix + "/part" +
                     std::to_string(active->rank() / num_per_part) + ".silo";
    pmpio.dirname = "/ensemble" + std::to_string(ens.label());

    pmpio.comm = active->split(active->rank() / num_per_part);
    fs::mpido(*active, [&]() { fs::create_directories(prefix); });
  }

  active->barrier();

  std::vector<particle::load_t> loads;
  loads.reserve(particles.size());
  for (auto sp : particles) loads.push_back(particles[sp].size());

  ens.reduce_to_chief(mpi::by::SUM, loads.data(), loads.size());

  pmpio([&](auto& dbfile) {
    // write global data
    if (!dbfile.exists("/timestep")) {
      dbfile.write("/cartesian_topology",
                   reinterpret_cast<const int*>(ens.cart_topos.begin()),
                   {DGrid});

      using T = std::underlying_type_t<particle::species>;
      std::vector<T> sps;
      for (auto sp : particles) sps.push_back(static_cast<T>(sp));
      dbfile.write("/species", sps);

      dbfile.write("/timestep", timestep);
    }

    {  // write ensemble-wise data
      if (!dbfile.exists("label")) {
        dbfile.write("cartesian_coordinates", ens.cart_coords.begin(), {DGrid});
        dbfile.write("label", ens.label());
      }

      FieldCkpt<R, DGrid> ckpt;
      if (ens.intra.rank() == ens.chief) {
        int idx = 0;
        for (auto sp : particles)
          dbfile.write(properties[sp].name + "_load", loads[idx++]);

        ckpt.save(dbfile, "E", E);
        ckpt.save(dbfile, "B", B);
      }
    }

    dbfile.mkcd("rank" + std::to_string(ens.intra.rank()));
    {  // write process specific data
      int r = ens.intra.rank();
      dbfile.write("r", r);
      ParticleArrayCkpt<R, S> ckpt;
      for (auto sp : particles) {
        ckpt.save(dbfile, properties[sp].name, particles[sp]);
      }
    }
    dbfile.cd("..");
  });
  active->barrier();  // ensure the checkpoint is fully saved before any other
                      // actions

  return prefix;
}

template <int DGrid, typename R, template <typename> class S>
int load_checkpoint(std::string dir,
                    std::optional<dye::Ensemble<DGrid>>& ens_opt,
                    const std::optional<mpi::CartComm>& cart_opt,
                    field::Field<R, 3, DGrid>& E, field::Field<R, 3, DGrid>& B,
                    particle::map<particle::array<R, S>>& particles,
                    const particle::map<particle::Properties>& properties,
                    int target_load) {
  int checkpoint_ts = 0;
  using T = std::underlying_type_t<particle::species>;
  int num_sps = 0;
  std::vector<T> sps;
  int num_ens = 0;
  particle::load_t myload = 0;

  // rank0 read data
  if (mpi::world.rank() == 0) {
    auto dir_itr = fs::directory_iterator(dir);
    auto dbfile = silo::open(*dir_itr, silo::Mode::Read);
    assert(dbfile.var_exists("/timestep"));
    assert(dbfile.var_exists("/cartesian_topology"));
    assert(dbfile.var_exists("/species"));

    dbfile.read("/timestep", &checkpoint_ts);
#if PIC_DEBUG
    lgr::file << "LDCKPT timestep = " << checkpoint_ts << std::endl;
    lgr::file << silo::errmsg() << std::endl;
#endif
    {
      int ndims = dbfile.var_length("/cartesian_topology");
      std::vector<int> topos(ndims);
      dbfile.read("/cartesian_topology", topos.data());
#if PIC_DEBUG
      lgr::file << "LDCKPT cartesian_topology = (";
      for (auto x : topos) lgr::file << x << ",";
      lgr::file << ")" << std::endl;
      lgr::file << silo::errmsg() << std::endl;
#endif
      num_ens = 1;
      for (auto x : topos) num_ens *= std::abs(x);
    }
    {
      num_sps = dbfile.var_length("/species");
      sps.resize(num_sps);
      dbfile.read("/species", sps.data());
#if PIC_DEBUG
      lgr::file << "LDCKPT species = (";
      for (auto x : sps) lgr::file << x << ",";
      lgr::file << ")" << std::endl;
      lgr::file << silo::errmsg() << std::endl;
#endif
    }
  }
  // rank0 broadcast data
  {
    int buf[3] = {checkpoint_ts, num_sps, num_ens};
    mpi::world.broadcast(0, buf, 3);
    if (mpi::world.rank() != 0) {
      checkpoint_ts = buf[0];
      num_sps = buf[1];
      num_ens = buf[2];
      sps.resize(num_sps);
    }
    mpi::world.broadcast(0, sps.data(), sps.size());
#if PIC_DEBUG
    if (mpi::world.rank() != 0) {
      lgr::file << "LDCKPT timestep = " << checkpoint_ts << std::endl;

      lgr::file << "LDCKPT num_ens = " << num_ens << std::endl;

      lgr::file << "LDCKPT species = (";
      for (auto x : sps) lgr::file << x << ",";
      lgr::file << ")" << std::endl;
    }
#endif
  }
  // each primary figure out its load
  if (cart_opt) {
    assert(ens_opt);
    const int mylabel = ens_opt->label();
    for (auto f : fs::directory_iterator(dir)) {
      bool is_found = false;
      auto sf = silo::open(f, silo::Mode::Read);
      for (const auto& dname : sf.toc_dir()) {
        if (sf.var_exists(dname + "/rank0")) {
          int l = sf.read1<int>(dname + "/label");
          if (mylabel == l) {
            sf.cd(dname);
            for (int i = 0; i < num_sps; ++i) {
              auto sp = static_cast<particle::species>(sps[i]);
              std::string entry = properties[sp].name + "_load";
              assert(sf.var_exists(entry));
              myload += sf.read1<particle::load_t>(entry);
            }
#if PIC_DEBUG
            lgr::file << "LDCKPT Label = " << mylabel << ", Load = " << myload
                      << std::endl;
            lgr::file << silo::errmsg() << std::endl;
#endif
            is_found = true;
            sf.cd("..");
            break;
          }
        }
      }
      if (is_found) break;
    }
  }

  // dynamic assign computational resources
  {
    auto ens_opt_new = dye::deploy(myload, ens_opt, cart_opt, target_load);
    ens_opt.swap(ens_opt_new);
  }

  // by this step, all processes should have ens_opt set properly. The real
  // loading begins now
  if (!ens_opt) return checkpoint_ts;

  const int mylabel = ens_opt->label();
  const int myrank = ens_opt->intra.rank();
  const auto ens_size = ens_opt->intra.size();

  FieldCkpt<R, DGrid> f_ckpt;
  ParticleArrayCkpt<R, S> p_ckpt;
  for (auto f : fs::directory_iterator(dir)) {
#if PIC_DEBUG
    lgr::file << "Reading file " << f << std::endl;
#endif
    auto sf = silo::open(f, silo::Mode::Read);
    for (const auto& dname : sf.toc_dir()) {
      if (dname.find("ensemble") != 0) continue;
      int l = sf.read1<int>(dname + "/label");
      if (l != mylabel) continue;
#if PIC_DEBUG
      lgr::file << "Reading EB from " << dname << std::endl;
#endif
      sf.cd(dname);
      if (0 == myrank &&
          sf.var_exists(
              "rank0")) {  // NOTE: ranks of one ensemble may exist across files
        f_ckpt.load(sf, "E", E);
        f_ckpt.load(sf, "B", B);
      }
      for (const auto& rdir : sf.toc_dir()) {
        if (rdir.find("rank") != 0) continue;
        sf.cd(rdir);
        int r = sf.read1<int>("r");
#if PIC_DEBUG
        lgr::file << "LDCKPT rank = " << r << std::endl;
        lgr::file << silo::errmsg() << std::endl;
#endif
        for (auto i : sps) {
          auto sp = static_cast<particle::species>(i);
          p_ckpt.load(sf, properties[sp].name, particles[sp], myrank + r,
                      ens_size);
        }
        sf.cd("..");
      }
      sf.cd("..");
    }
  }

  return checkpoint_ts;
}
}  // namespace ckpt
