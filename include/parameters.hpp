#ifndef  _PARAMETERS_H_
#define  _PARAMETERS_H_

#include <array>
#include <string>

// TODOL traits and params are basically the same thing. Just compile time known or not
template < typename Real, std::size_t DGrid >
struct Params {
  Real dt;
  int total_timesteps;

  Real e; // electric charge

  std::string this_run_dir;

  // Grid grid; // local simulation grid
  // Grid dataGrid; // local data export grid

  // std::array<bool, 3> is_periodic; // only needed in sendcellsleftright

  // the following are ensemble specs, which will be stored on primary and be passed on to all replicas
  struct Locale {
    int label;
    std::array< int, DGrid > anchor;
    std::array< int, DGrid > extent;
    std::array< std::array<Real, 2>, DGrid > borders;

    std::array< int, DGrid > coordinate;
    std::array< std::array<bool, 2>, DGrid > is_at_boundary;
    std::array< std::array<bool, 2>, DGrid > is_axis;
    std::array< std::array<int, 2>, DGrid > neighbors;
  } locale;

  // void Init(const Dashboard& db, const MPICommunicator& comm);

  // void SyncOverEnsemble(const MPIEnsembleCommunicator& ensemble);
};

#endif
