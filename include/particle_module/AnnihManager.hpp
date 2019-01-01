#ifndef _ANNIHMANAGER_H_
#define _ANNIHMANAGER_H_

#include "Particles.h"
#include "Fields.h"
#include "Dashboard.h"

class MPIEnsembleCommunicator;

struct AnnihManager {
  void MarkPairs( Particles<QCarrier>& electrons, Particles<QCarrier>& positrons, Scalar dt, const Grid& grid, const MPIEnsembleCommunicator& ensemble, const DBPane_Annihilation::Matter& pane );

  void Annihilate( Particles<QCarrier>& electrons, Particles<QCarrier>& positrons );
  void Annihilate( Particles<QNeutral>& photons, const Grid& grid, const DBPane_Annihilation::Photon& pane );
};

#endif
