#include "particle_module/ParticlePusher.hpp"
#include "core/particle.hpp"
#include "core/field.hpp"
#include "core/coordinate.hpp"

namespace particle :: force {
  void rad_cooling(Vec3<double>& p, const Vec3<double>& B, const Vec3<double>& E, double re, double dt) {
    double gamma = sqrt(1.0 + p.dot(p));
    auto beta = p / gamma;
    // std::cout << "beta is " << beta << std::endl;
    auto tmp = E + beta.cross(B);
    // std::cout << "tmp is " << tmp << std::endl;
    auto dp = tmp.cross(B) + E * beta.dot(E);
    // std::cout << "-------" << std::endl;
    // std::cout << "first term is " << dp << std::endl;
    double bE = beta.dot(E);
    auto sec = beta * (gamma * gamma * (tmp.dot(tmp) - bE * bE));
    // std::cout << "second term is " << sec << std::endl;
    dp -= sec;
    dp *= (2.0 * re / 3.0) * dt;
    // std::cout << "dp is " << dp << std::endl;
    p += dp;
  }

  void sync_cooling(Vec3<double>& p, const Vec3<double>& B, const Vec3<double>& E, double re, double dt) {
    // double gamma = sqrt(1.0 + p.dot(p));
    double B2 = B.dot(B);
    auto p_perp = p - B * (p.dot(B)) / B2;
    double pp = sqrt(p_perp.dot(p_perp));
    // auto dp
    // std::cout << "dp is " << dp << std::endl;
    p -= p_perp * (2.0 * re * dt * pp * B2 / 3.0);
    // p = B * (p.dot(B)) / B2;
  }

  void landau0_cooling(Vec3<double>& p, const Vec3<double>& B, const Vec3<double>& E) {
    Scalar EB2 = E.dot(B);
    EB2 = EB2 * EB2;
    Scalar B2_E2 = B.dot(B) - E.dot(E);
    // calculate E'^2
    Scalar Ep2 = 2 * EB2 / ( std::sqrt(B2_E2 * B2_E2 + 4 * EB2) + B2_E2 );
    auto beta_ExB = E.cross(B) / (B.dot(B) + Ep2);
    // find B' modulo gamma_ExB
    auto Bp = B - beta_ExB.cross(E);
    // obtain the momentum with perpendicular components damped
    p = Bp * ( p.dot(Bp) / Bp.dot(Bp) );
    p += beta_ExB * std::sqrt( ( 1.0 + p.dot(p) ) / ( 1.0 - beta_ExB.dot(beta_ExB) ) );
  }

  inline Vec3<MOM_TYPE> lorentz( Scalar dt, const Vec3<MOM_TYPE>& p, const Scalar& charge_to_mass, const Vec3<Scalar>& EVector, const Vec3<Scalar>& BVector ) {
    Scalar lambda =  charge_to_mass * dt / 2.0; // measured in units of (e/m) * R_* / c;
    Vec3<Scalar> u_halfstep = p + EVector * lambda + p.cross( BVector ) * ( lambda / std::sqrt( 1.0 + p.dot(p) ) );
    Vec3<Scalar> upr = u_halfstep + EVector * lambda;
    Vec3<Scalar> tau = BVector * lambda;
    // store some repeatedly used intermediate results
    Scalar tt = tau.dot( tau );
    Scalar ut = upr.dot( tau );

    Scalar sigma = 1.0 + upr.dot(upr) - tt;
    // inv_gamma2 means ( 1 / gamma^(i+1) ) ^2
    Scalar inv_gamma2 =  2.0 / ( sigma + std::sqrt( sigma * sigma + 4.0 * ( tt + ut * ut ) ) );
    Scalar s = 1.0 / ( 1.0 + inv_gamma2 * tt );
    Vec3<MOM_TYPE> p_vay = ( upr + tau * ( ut * inv_gamma2 ) + upr.cross( tau ) * std::sqrt(inv_gamma2) ) * s;
    //    Vec3<MOM_TYPE> dp = p_vay - p;
    return p_vay - p;
  }

  float CalculateRc( Scalar dt, const Vec3<Real> &p, const Vec3<Real> &dp ) {
    // find momentum at half time step
    Vec3<MOM_TYPE> phalf( p + dp * 0.5 );
    Vec3<MOM_TYPE> v( phalf / std::sqrt( 1.0 + phalf.dot(phalf) ) );
    Scalar vv = v.dot( v );
    Vec3<MOM_TYPE> a( dp / dt ); // a is for now force, will be converted to dv/dt
    // convert a to dv/dt
    a = ( a - v * ( v.dot(a) ) ) * std::sqrt( 1.0 - vv );
    Scalar va = v.dot( a ); // get the real v dot a
    return vv / std::max( std::sqrt( a.dot(a) - va * va / vv ), 1e-6 ); // in case denominator becomes zero
  }

  float GetDipolarRc( const Scalar &r_sph, const Scalar &cos_th, const Scalar &phi) {
    Scalar sin_th = std::sqrt( 1.0 - cos_th * cos_th );
    Scalar tmp1 = 1.0 + cos_th * cos_th;
    Scalar tmp2 = 3.0 * tmp1 - 2.0;
    return r_sph * tmp2 * std::sqrt(tmp2) / ( 3.0 * tmp1 * sin_th );
  }

}

// TODO
namespace particle {
  void update_p( Particle& ptc, Real dt, const Vec3<Real>& E, const Vec3<Real>& B ) {
    Vec3<Real> dp;
    auto& p = ptc.p();

    // Apply Lorentz force
    if ( _pane.lorentz_On && ( !check_bit( ptc.flag, ParticleFlag::ignore_em ) ) ) {
      dp += force::lorentz( dt, p, particles.Attributes().chargeToMass, E, B );
    }

    const auto& q = ptc.q();

    if ( _pane.gravity_On )
      dp += _pane.gravity( dt, q[0], q[1], q[2] );

    // TODO put Rc in a separate container
    // // FIXME don't use a uniform number for Rc
    // // ptc.Rc = CalculateRc( dt, p, dp );
    // ptc.Rc = 1.0;

    p += dp;

    // FIXME add control in dashboard for rad_cooling
    // if ( particles.Attributes().isRadiate ) {
    //   rad_cooling(p, BVector, EVector, 2e-9, dt);
    // }

    // when B is strong enough, damp the perpendicular component of momentum
    // FIXME ions should also be affected right?
    // TODO should this go before p += dp
    if ( _pane.landau0_On && particles.Attributes().isRadiate ) {
      if ( B.dot(B) > _pane.B_landau0 * _pane.B_landau0 ) {
        landau0_cooling( p, B, E );
      }
    }
  }
}



namespace particle {
  template < class Coord >
  void push ( std::vector<Particle>& particles, Real dt, const Params& params,
              const VectorField<Real>& EField, const VectorField<Real>& BField ) {
    const auto& grid = params.grid;
    for ( auto& ptc : particles ) {
      if( ptc.IsEmpty(i) ) continue;

      // update q
      ptc.q() += dq;

      // TODO
      // handle_boundary( ptc is_at_boundary, is_axis );

      if ( !check_bit( ptc.flag, ParticleFlag::ignore_force ) ) {
        // TODO get E B
        // TODO different forces may need different filtering
        update_p( ptc, dt, E, B );
      }
      auto dx = Coord::geodesic_move( ptc.q(), ptc.p(), dt );

    }
  }

  template < class Coord >
  void push ( std::vector<Particle>& particles, Real dt, const Params& params ) {
    const auto& grid = params.grid;
    for ( auto& ptc : particles ) {
      if ( ptc.IsEmpty() ) continue;

      // need species
      auto dq = Coord::geodesic_move( ptc.q(), ptc.p(), dt );
      ptc.q() += dq;

      // TODO move this to pairProducer
      // ptc.path_left -= dt;

      // TODO
      // handle_boundary( ptc is_at_boundary, is_axis );
    }
  }

}
