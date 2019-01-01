#include "ParticleSorter.h"
#include <algorithm>
#include "MagicNumbers.h"
#include "Logger.h"
#include <chrono>
#include "Grid.h"

using namespace std;
using namespace std::chrono;

void InitializeIndex(vector<int>& index, int num) {
  for (int i = 0; i < num; ++i) {
    index[i] = i;
  }
}

ParticleSorter::ParticleSorter()
    : _numMax(0), _numTiles(0) {}

ParticleSorter::ParticleSorter(const Grid& grid, size_t numMaxPtc)
    : _numMax(numMaxPtc),
      _tile(_numMax),
      _index(_numMax) {
  _numTiles = std::max(grid.tileNum(), 27) + 1; // see header file for definition. +1 for empty particles.
  _numInTile = vector<int>( _numTiles  );
  _scanTile = vector<int>(_numTiles );
}

ParticleSorter::~ParticleSorter() {}


template< typename P >
void ParticleSorter::Sort(Particles<P>& particles, int num, const Grid& grid, bool isCommMode) {
    //  InitializeIndex(_index, num);

    unsigned int sorttilenum = isCommMode ? 27 : grid.tileNum();
    sorttilenum += 1; // for storage of empty particles
    // Zero out _numInTile
    for( unsigned int i=0; i < _numTiles; i++){
        _numInTile[i] = 0;
    }

    // Figure out the indices of each particle
    for (int i = 0; i < num; ++i) {
      int cell = particles.PtcData()[i].cell;
      int c1 = grid.getC1(cell);
      int c2 = grid.getC2(cell);
      int c3 = grid.getC3(cell);

      if ( particles.IsEmpty(i) ) {
        _tile[i] = sorttilenum - 1; // put empty particles in the last tile
      }
      else if(isCommMode){
          int d1 = (c1 >= grid.guard[0]) + (c1 >= (grid.dims[0] - grid.guard[0]));
          int d2 = (c2 >= grid.guard[1]) + (c2 >= (grid.dims[1] - grid.guard[1]));
          int d3 = (c3 >= grid.guard[2]) + (c3 >= (grid.dims[2] - grid.guard[2]));
          _tile[i] = d1 + d2 * 3 + d3 * 9;
      }
      else{
          _tile[i] = grid.tileId(c1,c2,c3);
      }

      _index[i] = _numInTile[_tile[i]];
      _numInTile[_tile[i]] += 1;
    }

    // need to initialize _scanTile[0] ? Yes, in the exclusive, we want each _scanTile[i] to give the
    // the subscript of the first element in each bin.
    // This also includes empty particles.
    _scanTile[0] = 0;
    for (unsigned int i = 1; i < sorttilenum; ++i) {
      _scanTile[i] = _scanTile[i - 1] + _numInTile[i - 1];
    }

    for (int i = 0; i < num; ++i) {
      _index[i] += _scanTile[_tile[i]];
    }

    particles.Rearrange( _index, num );

    //when using particles.number() as num,
    // and when there are empty particles involved, adjust num
    if( particles.Number() == num && _numInTile[sorttilenum-1] != 0 )
        particles.SetNum( _scanTile[sorttilenum-1] );
    Logger::print(0, "Sort complete, we now have", particles.Number(), particles.NameStr(), "in the pool!");
}


template void ParticleSorter::Sort<QCarrier>(Particles<QCarrier> &particles, int num, const Grid& grid, bool isCommMode);
template void ParticleSorter::Sort<QNeutral>(Particles<QNeutral> &photons, int num, const Grid& grid, bool isCommMode);
