#ifndef  _PARTICLEUPDATER_HPP_
#define  _PARTICLEUPDATER_HPP_

#include "dashboard.hpp"

struct Params;
namespace mpi { struct Comm; }
struct DynamVars;
class Rng;

class ParticlePusher;
class ParticleInjector;
class CurrentDepositer;
class ParticleSorter;
class ParticleCommunicator;
class PairProducer;
struct AnnihManager;

namespace particle {
  void update(DynamVars& data, int timestep, const Params& params, const mpi::Comm& comm, Rng& rng);

}


template <CoordType Coord>
class ParticleUpdater {
public:
  ParticleUpdater(const Dashboard& db, const Params& params, const MPICommunicator& comm, std::unordered_map<ParticleType, ParticleSorter*>& pSorter, ParticleCommunicator*& ptc_comm);
  ~ParticleUpdater();
  void Update(AperData& data, int timestep, const AperParams& params, const MPICommunicator& comm, Rng& rng, Domain* domain);

private:
  int _sort2ndInterval;
  DBPane_Annihilation _annih_pane;

  ParticlePusher* _pPusher = nullptr;
  CurrentDepositer* _pDepositer = nullptr;
  std::unordered_map<BoundaryPosition, ParticleInjector* > _pInjector;
  PairProducer* _pPairProducer = nullptr;
  AnnihManager* _pAnnih = nullptr;

  std::unordered_map<ParticleType, ParticleSorter*>& _pSorter;
  ParticleCommunicator*& _ptc_comm;
};

#endif
