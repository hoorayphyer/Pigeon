#include <fmt/format.h>

#include <algorithm>
#include <stdexcept>

#include "filesys/filesys.hpp"
#include "simulator/argparser.hpp"
#include "simulator/builder.hpp"

namespace {

std::string init_this_run_dir(std::string prefix, std::string dirname) {
  // use world root time to ensure uniqueness
  std::string this_run_dir;
  if (mpi::world.rank() == 0) {
    prefix = fs::absolute(prefix);
    fs::remove_slash(prefix);
    fs::remove_slash(dirname);

    // in case of running too frequently within a minute, directories with
    // postfixed numbers are created
    if (fs::exists(prefix + "/" + dirname)) {
      for (int n = 1;; ++n) {
        if (!fs::exists(prefix + "/" + dirname + "-" + std::to_string(n))) {
          dirname += "-" + std::to_string(n);
          break;
        }
      }
    }
    this_run_dir = prefix + "/" + dirname;

    // local directory for storing data symlinks
    // TODO Maybe deprecate this?
#ifdef APPARENT_DATA_DIR
    std::string local_data_dir = []() {
      std::string str = APPARENT_DATA_DIR;
      fs::remove_slash(str);
      return str;
    }();
#else
    std::string local_data_dir = "Data";
#endif

    fs::create_directories(this_run_dir);
    fs::create_directories(local_data_dir);
    fs::create_directory_symlink(this_run_dir, local_data_dir + "/" + dirname);
  }

  if (mpi::world.size() > 1) {
    char buf[200];
    if (mpi::world.rank() == 0) {
      for (int i = 0; i < this_run_dir.size(); ++i) buf[i] = this_run_dir[i];
      buf[this_run_dir.size()] = '\0';
      mpi::world.broadcast(0, buf, 200);
    } else {
      mpi::world.broadcast(0, buf, 200);
      this_run_dir = {buf};
    }
  }

  return this_run_dir;
}

std::string data_dirname(std::string project_name) {
  char subDir[100] = {};
  for (int i = 0; i < 100; ++i) subDir[i] = '\0';
  if (mpi::world.rank() == 0) {
    char myTime[100] = {};
    time_t rawtime;
    struct tm* timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(myTime, 100, "%Y%m%d-%H%M", timeinfo);
    snprintf(subDir, sizeof(subDir), "%s", myTime);
  }
  return project_name + "-" + subDir;
}
}  // namespace

namespace pic {

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
SimulationBuilder<DGrid, R, S, RJ, RD>::SimulationBuilder(CLIArgs args)
    : m_args(std::move(args)) {
  auto& dir = m_args.resume_dir;
  if (dir) {
    *dir = fs::absolute(*dir);
  }

  mpi::initialize(m_args.rest);
}

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
SimulationBuilder<DGrid, R, S, RJ, RD>::~SimulationBuilder() {
  lgr::file.close();
  // TODO double check if this is OK
  mpi::world.barrier();
  // reset m_sim to force destruction of potential mpi communicators before
  // mpi::finalize
  m_sim.reset();
  mpi::finalize();
}

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
SimulationBuilder<DGrid, R, S, RJ, RD>&
SimulationBuilder<DGrid, R, S, RJ, RD>::initialize_this_run_dir(
    std::string prefix, std::string project_name) {
  auto dirname = data_dirname(project_name);
  auto this_run_dir = init_this_run_dir(prefix, dirname);

  m_this_run_dir.emplace(this_run_dir);

  populate_this_run_dir();

  auto log_file =
      fmt::format("{}/logs/rank/{}.log", *m_this_run_dir, mpi::world.rank());
  lgr::file.set_filename(log_file);

  return *this;
}

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
SimulationBuilder<DGrid, R, S, RJ, RD>&
SimulationBuilder<DGrid, R, S, RJ, RD>::create_cartesian_topology(
    const apt::array<int, DGrid>& dims,
    const apt::array<bool, DGrid>& periodic) {
  m_cart.emplace();
  const int nprmy =
      std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<int>());

  if (nprmy > mpi::world.size())
    throw std::invalid_argument(
        "Size of the Cartesian topology exceeds the size of world!");

  // simply use first few members in mpi::world as primaries
  bool is_prmy = mpi::world.rank() < nprmy;
  auto comm_tmp = mpi::world.split({is_prmy}, mpi::world.rank());
  if (is_prmy) {
    std::vector<int> d(DGrid);
    std::vector<bool> p(DGrid);
    for (int i = 0; i < DGrid; ++i) {
      d[i] = dims[i];
      p[i] = periodic[i];
    }
    m_cart->emplace(*comm_tmp, d, p);
  }

  m_dims = dims;
  m_periodic = periodic;

  return *this;
}

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
std::string SimulationBuilder<DGrid, R, S, RJ, RD>::precondition() const {
  std::string msg = "";
  if (!m_supergrid) {
    msg += "- must set supergrid\n";
  }
  if (!m_cart) {
    msg += "- must initialize cartesian topology\n";
  }
  if (!m_this_run_dir) {
    msg += "- must set this_run_dir\n";
  }
  if (m_ic_actions.size() == 0) {
    msg += "- no initial condition is given\n";
  }
  if (!m_fld_guard) {
    msg += "- must set the field guard\n";
  }
  return msg;
}

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
int SimulationBuilder<DGrid, R, S, RJ, RD>::load_init_cond(Simulator_t& sim) {
  if (mpi::world.rank() == 0)
    std::cout << "Loading initial condition..." << std::endl;
  int init_ts = 0;
  // TODO profiling_plan
  // if (profiling_plan.on and profiling_plan.is_qualified()) {
  //   lgr::file << "--- Loading Initial Condition..." << std::endl;
  // }

  if (m_args.resume_dir) {
    init_ts = ckpt::load_checkpoint(*m_args.resume_dir, sim.m_ens_opt,
                                    sim.m_cart_opt, sim.m_E, sim.m_B,
                                    sim.m_particles, sim.m_properties);
    if (sim.m_ens_opt) {
      sim.update_parts();

      auto bundle = typename PostResumeAction_t::Bundle_t{
          sim.m_E,    sim.m_B,       sim.m_J, sim.m_particles, sim.m_properties,
          sim.m_grid, sim.m_ens_opt, init_ts, *m_this_run_dir};

      for (auto& act : m_post_resume_actions) {
        // TODO also check range emptiness
        if (!act) continue;
        sim.taylor(act->ranges());
        (*act)(bundle);
      }
    }
    ++init_ts;  // checkpoint is saved at the end of a timestep
  } else {
    // FIXME a temporary fix, which may crash under the edge case in which
    // initially many particles are created.
    if (sim.m_cart_opt) {
      auto bundle = typename InitialConditionAction_t::Bundle_t{
          sim.m_E,         sim.m_B,          sim.m_J,
          sim.m_particles, sim.m_properties, sim.m_grid};
      for (auto& act : m_ic_actions) {
        // TODO also check range emptiness
        if (!act) continue;
        sim.taylor(act->ranges());
        (*act)(bundle);
      }

      for (auto sp : sim.m_particles) {
        for (auto ptc : sim.m_particles[sp]) {  // TODOL semantics
          ptc.set(sp);
        }
      }
      field::copy_sync_guard_cells(sim.m_E, *sim.m_cart_opt);
      field::copy_sync_guard_cells(sim.m_B, *sim.m_cart_opt);
      field::copy_sync_guard_cells(sim.m_J, *sim.m_cart_opt);
    }
  }
  // broadcast to ensure uniformity
  mpi::world.broadcast(0, &init_ts, 1);

  // TODO the following can be made into a separate action
  // // initialize trace_counters to ensure unique trace serial numbers across
  // // runs
  // pic::Traman::init(mpi::world.rank(), _particles);

  // TODO profiling plan
  // if (profiling_plan.on and profiling_plan.is_qualified()) {
  //   lgr::file << "--- Done." << std::endl;
  // }

  return init_ts;
}

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
void SimulationBuilder<DGrid, R, S, RJ, RD>::set_up_journal() const {
  const auto& this_run_dir = *m_this_run_dir;
  const std::string jnl_default(this_run_dir + "/journal.txt");
  std::string jnl = [&jnl_default, this]() {
    std::string res = jnl_default;
    if (m_args.journal_file) {
      // a journal file is specified
      res = fs::absolute(*m_args.journal_file);
      if (!fs::exists(res)) {
        // TODO std cout?
        std::cout << "Specified journal doesn't exist. Using default journal "
                     "instead."
                  << std::endl;
        res = jnl_default;
      }
    }

    return res;
  }();
  std::ofstream out;
  out.open(jnl, std::ios_base::app);  // NOTE app creates new file when jnl
                                      // doesn't exist
#if PIC_DEBUG
  out << "BuildType := Debug" << std::endl;
#else
  out << "BuildType := Release" << std::endl;
#endif
  out << "DataDir := " << this_run_dir << std::endl;
  if (m_args.resume_dir) out << "Resume := " << *m_args.resume_dir << std::endl;
  out.close();
  if (!fs::equivalent(jnl, jnl_default)) {
    // NOTE fs::rename doesn't work on some platforms because of
    // cross-device link.
    fs::copy_file(jnl, jnl_default);
  }
}

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
void SimulationBuilder<DGrid, R, S, RJ, RD>::populate_this_run_dir() const {
  fs::mpido(mpi::world, [this]() {
    const auto& this_run_dir = *m_this_run_dir;

    set_up_journal();

    fs::create_directories(this_run_dir + "/data");
    fs::create_directories(this_run_dir + "/logs");
    fs::create_directories(this_run_dir + "/pigeon");
    fs::copy_file("CMakeLists.txt", this_run_dir + "/pigeon/CMakeLists.txt");
    fs::copy_file("pic.hpp", this_run_dir + "/pigeon/pic.hpp");
    fs::copy_file("pic_impl.hpp", this_run_dir + "/pigeon/pic_impl.hpp");
    if (m_args.config_file) {
      fs::copy_file(*m_args.config_file, this_run_dir + "/pigeon/conf.toml");
    }
  });
}

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
SimulationBuilder<DGrid, R, S, RJ, RD>::Simulator_t&
SimulationBuilder<DGrid, R, S, RJ, RD>::build() {
  if (m_sim) {
    return *m_sim;
  }

  {
    auto msg = precondition();
    if (!msg.empty()) {
      msg = "Error: incomplete simulation setup\n" + msg;
      throw std::runtime_error(msg);
    }
  }

  // Simulator_t ctor is not accessible in optional
  m_sim.emplace(Simulator_t());

  auto& sim = *m_sim;
  sim.initialize(std::move(*m_supergrid), std::move(*m_cart),
                 std::move(m_props), std::move(m_periodic), m_dims);
  auto init_ts = load_init_cond(sim);

  sim.m_rng.set_seed(init_ts + mpi::world.rank());

  // TODO field/ptc actions and their orig ranges
  sim.m_f_prior_export = std::move(m_f_prior_export);
  sim.m_exporters = std::move(m_exporters);
  sim.m_f_post_export = std::move(m_f_post_export);
  sim.m_f_custom_step = std::move(m_f_custom_step);
  sim.m_fld_guard = *m_fld_guard;
  sim.m_this_run_dir = *m_this_run_dir;
  sim.m_initial_timestep = init_ts;

  return sim;
}
}  // namespace pic
