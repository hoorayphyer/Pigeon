#ifndef  _PARAMETERS_H_
#define  _PARAMETERS_H_

#include "types.hpp"

// TODOL traits and params are basically the same thing. Just compile time known or not
struct Params {
  Real dt;
  int total_timesteps;

  Real e; // electric charge
  Real ion_mass_ratio;
  Real ion_charge_ratio;

  std::string this_run_dir;
  // WeightType weightType;

  // Grid grid; // local simulation grid
  // Grid dataGrid; // local data export grid

  // std::array<bool, 3> is_periodic; // only needed in sendcellsleftright

  // the following are ensemble specs, which will be stored on primary and be passed on to all replicas
  // struct EnsembleSpecs {
  //   int label;
  //   std::array<int, 3> coordinate;
  //   std::array<bool, 6> is_at_boundary;
  //   std::array<bool, 6> is_axis;
  //   std::array<int, 3> neighbor_left;
  //   std::array<int, 3> neighbor_right;
  // } ens_specs;

  // void Init(const Dashboard& db, const MPICommunicator& comm);

  // void SyncOverEnsemble(const MPIEnsembleCommunicator& ensemble);
};

#endif
