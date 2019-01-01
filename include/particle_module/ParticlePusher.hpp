#ifndef _PARTICLEPUSHER_HPP_
#define _PARTICLEPUSHER_HPP_

template <typename T>
struct VectorField<T>;

// TODO check missing 4pi's maybe
namespace particle {

  // TODOL add species check for correct version of push
  template < class Coord >
  void push ( std::vector<Particle>& particles, Real dt, const Params& params,
              const VectorField<Real>& EField, const VectorField<Real>& BField );

  template < class Coord >
  void push ( std::vector<Particle>& particles, Real dt, const Params& params );
}



// #include "Particles.h"
// #include "Fields.h"
// #include "Dashboard.h"

// enum class CoordType;
// struct AperParams;

// class ParticlePusher {
// public:
//   ParticlePusher( const DBPane_Pusher& pane );

//   template < CoordType Coord >
//   void Push ( Particles<QCarrier>& particles, Scalar dt,
//               const VectorField<Scalar>& EField, const VectorField<Scalar>& BField,
//               const AperParams& params ) const;

//   template < CoordType Coord >
//   void Push ( Particles<QNeutral>& photons, Scalar dt, const AperParams& params ) const;

// private:
//   const DBPane_Pusher& _pane;

//   void update_x( Vec3<POS_TYPE>& x, Vec3<POS_TYPE>& dx, int& cell ) const;

//   void update_p ( Particles<QCarrier>& particles, Scalar dt, const Grid& grid,
//                        const VectorField<Scalar>& EField, const VectorField<Scalar>& BField,
//                        WeightType weight_type ) const;

//   template < CoordType Coord >
//   void update_dx (Particles<QCarrier>& particles, Scalar dt, const Grid& grid) const;

//   template < typename P >
//   void HandleBoundaries(Particles<P> &particles, std::array<bool, 6> is_bdry, std::array<bool, 6> is_axis) const;

//   // FIXME TODO: move Rc calculation to PairProducer?
//   float CalculateRc( Scalar dt, const Vec3<MOM_TYPE>& p, const Vec3<MOM_TYPE>& dp ) const;

//   // FIXME TODO: too specific to pulsar, temporarily defaulted to Spherical implementation.
//   float GetDipolarRc( const Scalar& x, const Scalar& y, const Scalar& z ) const;
// };

#endif
