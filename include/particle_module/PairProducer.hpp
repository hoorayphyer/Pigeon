#ifndef _PAIRPRODUCER_H_
#define _PAIRPRODUCER_H_

#include "Types.h"
#include "Dashboard.h"

class AperData;
class Rng;
struct Grid;

class PairProducer {
private:
  virtual void ProducePairsAndPhotonsImpl( int timestep, AperData& data, Scalar dt, const Grid &grid, Rng& rng  ) const {}

public:
  PairProducer() {}
  virtual ~PairProducer() {}

  inline void ProducePairsAndPhotons( int timestep, AperData& data, Scalar dt, const Grid &grid, Rng& rng ) const {
    ProducePairsAndPhotonsImpl( timestep, data, dt, grid, rng );
  }

}; // ----- end of class PairProducer

class PairProducer_Instant : public PairProducer {
private:
  void ProducePairs( int timestep, AperData& data, Scalar dt, const Grid &grid, Rng& rng ) const;

  virtual void ProducePairsAndPhotonsImpl( int timestep, AperData& data, Scalar dt, const Grid &grid, Rng& rng ) const override {
    ProducePairs( timestep, data, dt, grid, rng );
  }

protected:
  DBPane_PairProduction _pane;

  inline Scalar SamplePhotonEnergy( const Scalar& gamma, const float& Rc ) const;
  inline bool IsProduceCurvPhoton(Scalar dt, const Scalar& gamma, const float& Rc, Scalar randnum ) const;

public:
  PairProducer_Instant( const DBPane_PairProduction& pane );
  virtual ~PairProducer_Instant() override {}

}; // ----- end of class PairProducer_Instant : public PairProducer

class PairProducer_Photon : public PairProducer_Instant {
private:
  inline bool IsMagneticConvert( Scalar dt, const Vec3<Scalar>& Bvec, Scalar randnum ) const; // check the condition for photons to magnetic convert

  void ProducePairs( int timestep, AperData& data, Scalar dt, const Grid& grid, Rng& rng ) const;
  void ProducePhotons( int timestep, AperData& data, Scalar dt, const Grid &grid, Rng& rng ) const;
  // for photon channels, there will be two steps.
  void ProducePairsAndPhotonsImpl( int timestep, AperData& data, Scalar dt, const Grid &grid, Rng& rng) const override {
    ProducePairs( timestep, data, dt, grid, rng );
    ProducePhotons( timestep, data, dt, grid, rng );
  }

public:
  PairProducer_Photon( const DBPane_PairProduction& pane );
  ~PairProducer_Photon() override {}

}; // ----- end of class PairProducer_Photon : public PairProducer_Instant

#endif // ----- end of _PAIRPRODUCER_H_
