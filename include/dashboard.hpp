#ifndef _DASHBOARD_HPP_
#define _DASHBOARD_HPP_

// #include "Types.h"
// #include <array>
// #include <vector>
// #include <unordered_set>
// #include <unordered_map>
// #include <memory>
// #include <chrono>
// #include "CoordSystem.h"

// class InitialCondition;

// struct DBPane_Spacetime {
//   CoordType coordType;
//   int num_spatial_dimensions;
//   Scalar time_span;

//   Scalar q1_min;
//   Scalar q1_max;
//   bool periodic1;

//   Scalar q2_min;
//   Scalar q2_max;
//   bool periodic2;

//   Scalar q3_min;
//   Scalar q3_max;
//   bool periodic3;
// };

// struct DBPane_Particles {
//   Scalar e; // unit electric charge
//   WeightType weightType;

//   std::unordered_map<ParticleType, Scalar> mass;
//   std::unordered_map<ParticleType,size_t> capacity;

//   auto species() const {
//     std::unordered_set<ParticleType> sp;
//     for ( const auto& elm : capacity )
//       sp.emplace( elm.first );
//     return sp;
//   }

// };

// struct DBPane_Discretization {
//   Scalar dt;

//   int superN1; // doesn't include guard cells
//   int superN2;
//   int superN3;

//   int guard1;
//   int guard2;
//   int guard3;

//   int sampleN1; // doesn't include guard cells
//   int sampleN2;
//   int sampleN3;
// };

// struct DBPane_BoundaryCondition {
//   struct FBC {
//     FieldBCType type;
//     int indent;

//     // FIXME how to better set boundaryCondition without exposing too much unneeded information

//     // for rotating conductor
//     BCFunc_split_t ft;
//     BCFunc_split_x B1;
//     BCFunc_split_x B2;
//     BCFunc_split_x B3;
//     BCFunc_split_x E1;
//     BCFunc_split_x E2;
//     BCFunc_split_x E3;

//     // for damping
//     Scalar damping_rate;
//     // Scalar thickness; // physical value FIXME setting this requires coordinate system
//     // std::function<Scalar(Scalar, Scalar)> profile;
//   };
//   struct PBC {
//     ParticleBCType type;

//     // for injection
//     Scalar N_inj;
//     ParticleType q_minus;
//     ParticleType q_plus;
//     Profile_Func x_profile;
//     Momentum_XFunc p_profile;
//     Vec3<Scalar> sigma_p;
//     Scalar J_reg_x; // current regulated multiplicity
//     Scalar t_moving_average; // the physical time over which to average the current
//     Scalar layer;
//   };

//   std::unordered_map< BoundaryPosition, FBC > fieldBC;
//   std::unordered_map< BoundaryPosition, PBC > ptcBC;
//   std::unordered_set< BoundaryPosition > bdry_axissymmetry;
// };

struct DBPane_DataExport {
  int interval;
  int numPmpio;
};

// struct DBPane_Snapshot {
//   using Clock = std::chrono::high_resolution_clock;
// #define DUR(x) std::chrono::duration_cast<Clock::duration>(std::chrono::hours(x))

//   bool isOn;
//   int timestep_interval;
//   int numSnapshots;
//   int numPmpio;
//   int wallclock_interval; // in wall clock hours, setting value < 1 turns it off

//   bool isSave( int timestep ) const {
//     if ( !isOn ) return false;
//     if ( 0 == timestep ) return false;
//     if ( timestep % timestep_interval == 0 ) return true;
//     if ( wallclock_interval < 1 ) return false;
//     static auto t_last = Clock::now();
//     auto t = Clock::now();
//     if ( t - t_last > DUR(wallclock_interval) ) {
//       t_last = t;
//       return true;
//     }
//     return false;
//   }
// #undef DUR
// };

// struct DBPane_Parallelization {
//   std::array<int,3> partition;

//   struct Dynamic {
//     bool isOn;
//     int interval;
//     size_t targetLoadPerProc;
//   } dynamic;

//   std::vector< std::array<int,4> > replicaMap;
// };

// struct DBPane_Tracking {
//   bool isOn;
//   int startTimestep;
//   int endTimestep; // tracking new particles is stopped after endTimestep, but writing tracked particle data is not afftected
//   Scalar ratio;
//   int interval; // tracking interval needs to divide snapshot interval if restart safety is desired.
//   int numPmpio;
//   bool isReset; // whether to ignore any tracking information and start tracking all over
// };

// struct DBPane_Pusher {
//   bool lorentz_On;

//   bool landau0_On;
//   Scalar B_landau0;

//   bool gravity_On;
//   std::function<Vec3<MOM_TYPE>(Scalar, Scalar, Scalar, Scalar)> gravity;
// };

// struct DBPane_PairProduction {
//   PairProductionScheme pairScheme;
//   Scalar gamma_off; // energy under which emission is prohibited.
//   Scalar K_thr; // prefactor in calculating Gamma_thr
//   Scalar E_ph; // energy of curvature photon in units of mc^2
//   Scalar K_curv_em_rate; // prefactor for calculating emission rate

//   Scalar l_coll; // free path for gamma ray photon to collide with a X ray photon
//   Scalar B_magconv; // the threshold of magnetic field beyond which magnetic conversion is turned on
//   Scalar l_magconv; // effective free path for photons undergoing magnetic conversion.

//   std::function<bool( Scalar, Scalar, Scalar)> is_ignore_pair_production;
// };

// struct DBPane_Annihilation {
//   struct Matter{
//     bool isOn;
//     int interval;
//     int num_pairs_target;
//     Scalar rate;
//     bool (*annihilable) (Scalar, Scalar, Scalar);
//   } matter;

//   struct Photon{
//     bool (*annihilable) (Scalar, Scalar, Scalar);
//   } photon;
// };

// struct Dashboard {
//   std::string storageDir;
//   int sort2ndInterval;

//   DBPane_Spacetime spacetime;
//   DBPane_Particles particles;
//   DBPane_Discretization discretization;
//   DBPane_BoundaryCondition boundaryCondition;

//   DBPane_Pusher pusher;
//   DBPane_PairProduction pairProduce;
//   DBPane_DataExport dataExport;
//   DBPane_Snapshot snapshot;
//   DBPane_Annihilation annih;
//   DBPane_Parallelization parallel;
//   DBPane_Tracking track;

//   InitialCondition* initialCondition = nullptr;

//   ~Dashboard() { delete initialCondition; }
// };

#endif // DASHBOARD_H
