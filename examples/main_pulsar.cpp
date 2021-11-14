#include <cassert>
#include <iostream>

#include "pigeon.hpp"

// TODOL users may forget to sync value and state. Add another layer then
namespace particle {  // must have this namespace for now
template <typename T>
struct Specs {
  using value_type = T;
  static constexpr int Dim = 3;
  using state_type = long long;
  static_assert(8 * sizeof(state_type) >= 64);
};
}  // namespace particle

constexpr int DGrid = 2;
using real_t = float;
using real_j_t = float;
using real_export_t = float;
using ShapeF = particle::shapef_t<particle::shape::Cloud_In_Cell>;

using pgn = PIGEON<DGrid, real_t, particle::Specs, real_j_t, real_export_t>;

struct SampleFieldAction : public pgn::FieldAction_t {
  void operator()(const Bundle_t& bundle) const override {}
};

auto set_up_particle_properties() {
  particle::map<particle::Properties> properties;
  return properties;
}

constexpr real_t compile_time_const = 0;

int main(int argc, char** argv) {
  const auto args = pic::parse_args(argc, argv);
  assert(args.config_file);

  auto conf = pgn::ConfFile_t::load(*args.config_file);

  auto datadir_prefix = conf["datadir_prefix"].as_or<std::string>("../Data/");
  auto project_name = conf["project_name"].as_or<std::string>("Unnamed");
  auto dt = conf["dt"].as<real_t>();
  auto n_timesteps = conf["total_timesteps"].as_or<int>(100);
  auto gamma_fd = conf["pairs"]["gamma_fd"].as<real_t>();

  // TODO put these to toml after getting rid of apt::array
  apt::array<int, DGrid> dims = {1, 1};
  apt::array<bool, DGrid> periodic = {false, false};

  pgn::SimulationBuilder_t builder(args);

  builder.initialize_this_run_dir(datadir_prefix, project_name)
      .create_cartesian_topology(dims, periodic);

  // TODO
  auto properties = set_up_particle_properties();
  builder.set_particle_properties(properties);

  if (mpi::world.rank() == 0)
    std::cout << "Initializing simulator..." << std::endl;

  {
    auto& field_act = builder.add_field_action<SampleFieldAction>();
    { field_act.set_name("blahblah"); }
    // ....
  }

  auto& sim = builder.build();
  // sim.start();
  // sim.print_steps(); // prints steps in order

  // TODO having to have user call this is a bit error_prone
  const auto init_ts = sim.initial_timestep();
  if (mpi::world.rank() == 0) std::cout << "Launch" << std::endl;
  for (int ts = init_ts; ts < init_ts + n_timesteps; ++ts) {
    sim.evolve(ts, dt);
  }

  return 0;
}
