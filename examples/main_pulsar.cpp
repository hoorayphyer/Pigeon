#include <iostream>

#include "pigeon.hpp"

struct SampleFieldAction : public FieldAction {};

constexpr real_t compile_time_const = 0;

int main(int argc, char** argv) {
  const auto args = pic::parse_args(argc, argv);
  assert(args.config_file);

  auto conf = ConfFile::load(*args.config_file);

  // TODO can we absorb this into parse_args and use std::exit therein?
  if (args.is_dry_run) {
    int retcode = 0;
    std::cout << "Dry Run Checks :=" << std::endl;
#if PIC_DEBUG
    std::cout << "\tDebug" << std::endl;
#else
    std::cout << "\tRelease" << std::endl;
#endif
    // TODO proofread
    // std::cout << pic::proofread("\t") << std::endl;
    if (args.resume_dir) {
      auto resume_dir = fs::absolute(*args.resume_dir);
      if (!fs::exists(resume_dir)) {
        retcode = 1;
        std::cout << "ERROR : Invalid resume directory : " << resume_dir
                  << std::endl;
      } else if (resume_dir.find("checkpoints/timestep") == std::string::npos) {
        // if the directory exists but is not one of the checkpoints
        retcode = 1;
        std::cout << "ERROR : Invalid resume directory : " << resume_dir
                  << ". Specify which checkpoint!" << std::endl;
      } else {
        std::cout << "\tResume from : " << resume_dir << std::endl;
      }
    }

    return retcode;
  }

  auto dt = conf["dt"].as<real_t>();
  auto n_timesteps = conf["total_timesteps"].as_or<int>(100);
  auto gamma_fd = conf["pairs"]["gamma_fd"].as<real_t>();

  SimulationBuilder builder(args);

  // TODO
  mpi::commit(
      mpi::Datatype<particle::Particle<pic::real_t, particle::Specs>>{});

  builder.initialize_this_run_dir().create_cartesian_topology();

  // TODO
  auto properties = pic::set_up_particle_properties();
  builder.set_particle_properties(properties).;

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
    sim.evolve(ts, pic::dt);
  }

  return 0;
}
