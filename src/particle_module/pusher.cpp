#include "particle_module/pusher.hpp"
#include "core/particle.hpp"
#include "core/field.hpp"
#include "core/coordinate.hpp"

namespace vn = vec::numerical;

namespace particle :: force {
  // void rad_cooling(Vec3<double>& p, const Vec3<double>& B, const Vec3<double>& E, double re, double dt) {
  //   double gamma = sqrt(1.0 + p.dot(p));
  //   auto beta = p / gamma;
  //   // std::cout << "beta is " << beta << std::endl;
  //   auto tmp = E + beta.cross(B);
  //   // std::cout << "tmp is " << tmp << std::endl;
  //   auto dp = tmp.cross(B) + E * beta.dot(E);
  //   // std::cout << "-------" << std::endl;
  //   // std::cout << "first term is " << dp << std::endl;
  //   double bE = beta.dot(E);
  //   auto sec = beta * (gamma * gamma * (tmp.dot(tmp) - bE * bE));
  //   // std::cout << "second term is " << sec << std::endl;
  //   dp -= sec;
  //   dp *= (2.0 * re / 3.0) * dt;
  //   // std::cout << "dp is " << dp << std::endl;
  //   p += dp;
  // }

  // void sync_cooling(Vec3<double>& p, const Vec3<double>& B, const Vec3<double>& E, double re, double dt) {
  //   // double gamma = sqrt(1.0 + p.dot(p));
  //   double B2 = B.dot(B);
  //   auto p_perp = p - B * (p.dot(B)) / B2;
  //   double pp = sqrt(p_perp.dot(p_perp));
  //   // auto dp
  //   // std::cout << "dp is " << dp << std::endl;
  //   p -= p_perp * (2.0 * re * dt * pp * B2 / 3.0);
  //   // p = B * (p.dot(B)) / B2;
  // }

  auto landau0 =
    [](const auto& p, const auto& E, const auto& B) {
      Real EB2 = vn::E.dot(B);
      EB2 = EB2 * EB2;
      Real B2_E2 = vn::abs_sq(B) - vn::abs_sq(E);
      // calculate E'^2
      Real Ep2 = 2 * EB2 / ( std::sqrt(B2_E2 * B2_E2 + 4 * EB2) + B2_E2 );
      auto beta_ExB = vn::cross(E,B) / ( vn::abs_sq(B) + Ep2);
      // find B' modulo gamma_ExB
      auto Bp = B - vn::cross( beta_ExB, E);
      // obtain the momentum with perpendicular components damped
      auto p_new = Bp * ( vn::dot( p, Bp ) / vn::abs_sq(Bp) );
      p_new += beta_ExB * std::sqrt( ( 1.0 + vn::abs_sq(p_new) ) / ( 1.0 - vn::abs_sq(beta_ExB) ) );
      return p_new - p;
    }

  auto lorentz = // lambda = dt / mass * e/m NOTE this is actually rescaling Lorentz force
    []( Real lambda, const auto& p, const auto& E, const auto& B ) noexcept {
      // TODO optimize use of intermediate variables
      lambda /= 2.0;

      auto u_halfstep = p + E * lambda + vn::cross(p, B) * ( lambda / std::sqrt( 1.0 + vn::abs_sq(p) ) );
      auto upr = u_halfstep + E * lambda;
      auto tau = B * lambda;
      // store some repeatedly used intermediate results
      Real tt = vn::abs_sq(tau);
      Real ut = vn::dot(upr, tau);

      Real sigma = 1.0 + vn::abs_sq(upr) - tt;
      Real inv_gamma2 =  2.0 / ( sigma + std::sqrt( sigma * sigma + 4.0 * ( tt + ut * ut ) ) ); // inv_gamma2 means ( 1 / gamma^(i+1) ) ^2
      Real s = 1.0 / ( 1.0 + inv_gamma2 * tt );
      auto p_vay = ( upr + tau * ( ut * inv_gamma2 ) + vn::cross(upr, tau) * std::sqrt(inv_gamma2) ) * s;
      return p_vay - p;
    };

  // float CalculateRc( Scalar dt, const Vec3<Real> &p, const Vec3<Real> &dp ) {
  //   // find momentum at half time step
  //   Vec3<MOM_TYPE> phalf( p + dp * 0.5 );
  //   Vec3<MOM_TYPE> v( phalf / std::sqrt( 1.0 + phalf.dot(phalf) ) );
  //   Scalar vv = v.dot( v );
  //   Vec3<MOM_TYPE> a( dp / dt ); // a is for now force, will be converted to dv/dt
  //   // convert a to dv/dt
  //   a = ( a - v * ( v.dot(a) ) ) * std::sqrt( 1.0 - vv );
  //   Scalar va = v.dot( a ); // get the real v dot a
  //   return vv / std::max( std::sqrt( a.dot(a) - va * va / vv ), 1e-6 ); // in case denominator becomes zero
  // }

  // float GetDipolarRc( const Scalar &r_sph, const Scalar &cos_th, const Scalar &phi) {
  //   Scalar sin_th = std::sqrt( 1.0 - cos_th * cos_th );
  //   Scalar tmp1 = 1.0 + cos_th * cos_th;
  //   Scalar tmp2 = 3.0 * tmp1 - 2.0;
  //   return r_sph * tmp2 * std::sqrt(tmp2) / ( 3.0 * tmp1 * sin_th );
  // }

}

// TODO move check of forces on_off to somewhere else
namespace particle {
  template < typename Ptc, typename Vector >
  Vec<Real,Ptc::Dim> update_p( Ptc& ptc, const Species& sp, Real dt, const Vector& E, const Vector& B ) {
    Vec<Real, Ptc::Dim> dp;

    // Apply Lorentz force
    if ( _pane.lorentz_On && ( !ptc.is<flag::ignore_em>() ) ) {
      dp += force::lorentz( dt/sp.mass, ptc.p, E, B );
    }

    if ( _pane.gravity_On )
      dp += _pane.gravity( ptc.q ) * dt; // TODO gravity interface

    // TODO put Rc in a separate container
    // // FIXME don't use a uniform number for Rc
    // // ptc.Rc = CalculateRc( dt, p, dp );
    // ptc.Rc = 1.0;

    // FIXME add control in dashboard for rad_cooling
    // if ( sp.is_radiative ) {
    //   rad_cooling(p, BVector, EVector, 2e-9, dt);
    // }

    // when B is strong enough, damp the perpendicular component of momentum
    // FIXME ions should also be affected right?
    // TODO should this go before p += dp
    if ( _pane.landau0_On && sp.is_radiative ) {
      if ( vn::abs_sq(B) > _pane.B_landau0 * _pane.B_landau0 ) {
        dp += force::landau0( ptc.p, E, B );
      }
    }

    ptc.p += dp;

    return dp;
  }

  template < typename Ptc, CoordSys CS >
  Vec<Real,Ptc::Dim> update_q( Ptc& ptc, const Species& sp, Real dt ) {
    Real gamma = std::sqrt( (sp.mass > 0) + vn::abs_sq(ptc.p) );

    if constexpr ( CS == CoordSys::Cartesian ) {
      return coord<CS>::geodesic_move( ptc.q, ptc.p, dt / gamma );
    } else {
      auto dq = coord<CS>::geodesic_move( ptc.q, (ptc.p /= gamma), dt );
      ptc.p *= gamma;
      return dq;
    }

  }


}


namespace particle {
  // template < class Coord >
  // void push ( std::vector<Particle>& particles, Real dt, const Params& params,
  //             const VectorField<Real>& EField, const VectorField<Real>& BField ) {
  //   const auto& grid = params.grid;
  //   for ( auto& ptc : particles ) {
  //     if( ptc.IsEmpty() ) continue;

  //     // update q
  //     ptc.q() += dq;

  //     // TODO
  //     // handle_boundary( ptc is_at_boundary, is_axis );

  //     if ( !check_bit( ptc.flag, ParticleFlag::ignore_force ) ) {
  //       // TODO get E B
  //       // TODO different forces may need different filtering
  //       update_p( ptc, dt, E, B );
  //     }
  //     auto dx = Coord::geodesic_move( ptc.q(), ptc.p(), dt );

  //   }
  // }

  // template < class Coord >
  // void push ( std::vector<Particle>& particles, Real dt, const Params& params ) {
  //   const auto& grid = params.grid;
  //   for ( auto& ptc : particles ) {
  //     if ( ptc.IsEmpty() ) continue;

  //     // need species
  //     auto dq = Coord::geodesic_move( ptc.q(), ptc.p(), dt );
  //     ptc.q() += dq;

  //     // TODO move this to pairProducer
  //     // ptc.path_left -= dt;

  //     // TODO
  //     // handle_boundary( ptc is_at_boundary, is_axis );
  //   }
  // }

}
