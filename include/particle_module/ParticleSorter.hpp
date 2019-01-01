#ifndef  _PARTICLESORTER_H_
#define  _PARTICLESORTER_H_

#include <utility>
#include "Particles.h"

struct Grid;

class ParticleSorter
{
 private:
  size_t _numMax;  // _numMax is max of numMaxParticles and numMaxPhotons
  size_t _numTiles; // _numTiles refers to max total number of sort tiles used in one call to sort.

  std::vector<int> _tile;
  std::vector<int> _index;
  std::vector<int> _numInTile;
  std::vector<int> _scanTile;


 public:
  ParticleSorter();
  ParticleSorter(const Grid& grid, size_t numMaxPtc);
  virtual ~ParticleSorter();

  // there are two modes of sorts, communication mode and parallel mode.
  // Communication mode sorts over entire node WITHOUT subdividing the bulk
  // into various tiles, i.e. the bulk is in only one sorttile in this mode.
  // Parallel mode is usually used after communication, i.e. guard cells are free of particles and
  // all particles are in the bulk. In this mode, the bulk is subdivided into various tiles so that
  // after the sort, nearby particles are also nearby stored, which in parallel computation increases
  // speed.
  // The parameter isCommMode distinguishes these two situations.
  // P can stands for Particles or Photons
  template< typename P >
  void Sort(Particles<P>& particles, int num, const Grid& grid, bool isCommMode);

  // data accessors
  const std::vector<int>& getScanTile() const { return _scanTile; }
  const std::vector<int>& getNumInTile() const { return _numInTile; }
  const std::vector<int>& getIndex() const { return _index;}
  const std::vector<int>& getTile() const { return _tile;}
  size_t getNumTiles() const { return _numTiles;}
}; // ----- end of class ParticleSorter -----


#endif   // ----- #ifndef _PARTICLESORTER_H_  ----- 
