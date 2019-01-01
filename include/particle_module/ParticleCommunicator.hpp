#ifndef PARTICLECOMMUNICATOR_H
#define PARTICLECOMMUNICATOR_H

#include "Particles.h"
class ParticleSorter;
class EnsParticleGatherer;
class MPICommunicator;
class MPIEnsembleCommunicator;
struct AperParams;

#define COUNT_COMMTAGSPARTICLES 18  // This is the count of VALID directions
                                    // defined in CommTagsParticles.
                                    // FIXME: Are there too many buffers needed??

class ParticleCommunicator {
private:
  size_t _sizeParticleBuffer;
  std::array<Particles<QCarrier>, COUNT_COMMTAGSPARTICLES> _bufferParticles;
  std::array<Particles<QNeutral>, COUNT_COMMTAGSPARTICLES> _bufferPhotons;
  int _ens_rank_import = 0; // the rank in ensemble that will import all received particles as a result of communication among cartesian. This variable will change every time step
  EnsParticleGatherer* _ens_gatherer;

  // return type specific references to send buffers
  template < typename P >
  inline std::array<Particles<P>,COUNT_COMMTAGSPARTICLES>& GetSendBuffer() {}

  template< typename P >
  void SendRecvDirectional(Particles<P> &particles, const ParticleSorter& sorter, const MPICommunicator& comm, const Vec3<int>& sendVec, int direction, int neighbor_left, int neighbor_right);

  // copy newly arrived particles in the CENTER buffer to particle array
  template< typename P >
  void ImportRecvedParticles(Particles<P>& particles, const MPIEnsembleCommunicator& ensemble);

public:
  ParticleCommunicator( const MPICommunicator& comm );
  ~ParticleCommunicator();

  template< typename P >
  void SendRecvParticles( int timestep, Particles<P>& particles, const ParticleSorter& sorter, const MPICommunicator& comm, const AperParams& params );

  // with updated domain, ens gather buffers on primaries will be resized
  void Reset( int ens_size );

  friend class TestParticleCommnicator;

};

#endif // PARTICLECOMMUNICATOR_H
