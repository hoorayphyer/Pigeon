#ifndef _FIELD_CARTESIAN_SOLVER_IMPL_HPP_
#define _FIELD_CARTESIAN_SOLVER_IMPL_HPP_

#include "field/cartesian_solver/updater.hpp"
#include "field/sync.hpp"

namespace field {
  template < typename T >
  constexpr T diff( const T& f, int stride ) noexcept {
    return *(&f + stride) - f;
  }

  template < typename T >
  constexpr T diff2( const T& f, int stride ) noexcept {
    return *(&f + stride) + *(&f - stride) - 2.0 * f;
  }

  // NOTE there is no cross terms like DxDy(f), otherwise we need a whole new interface

}

namespace field {
  // TODOL C++20 has templated lambda, following classes can be put this inside operator()
  template < int DGrid, typename Real >
  struct Derivative {
    const apt::array<int,3>& s; // stride
    const mani::Grid<Real,DGrid>& g; // grid

    template < int I, typename T >
    constexpr T D( const T& f ) const noexcept {
      if constexpr ( I >= DGrid ) return 0.0;
      else return field::diff( f, s[I] ) / g[I].delta();
    }

    template < int I, typename T >
    constexpr T DD( const T& f ) const noexcept {
      if constexpr ( I >= DGrid ) return 0.0;
      else return field::diff2( f, s[I] ) / (g[I].delta() * g[I].delta());
    }

    template < typename T >
    constexpr T Lapl( const T& f ) const noexcept {
      return DD<0>(f) + DD<1>(f) + DD<2>(f);
    }

    template < int I, typename T >
    constexpr T Curl( const T& f_J, const T& f_K ) const noexcept {
      return D<(I+1)%3>(f_K) - D<(I+2)%3>( f_J );
    }

  };

  template < int I, int DGrid, typename Real, typename Real2 >
  inline void add_plus_curl( Field<Real,3,DGrid>& f, Real pre, const Field<Real2,3,DGrid>& f2, const Derivative<DGrid,Real>& d ) {
    constexpr int J = (I+1)%3;
    constexpr int K = (I+2)%3;
    for ( int i = 0; i < f[0].data().size(); ++i ) {
      f[I][i] += pre * d.template Curl<I>( f2[K][i], f2[J][i] );
    }
  }


  template < typename Real, int DGrid, typename RealJ >
  void CartesianUpdater<Real,DGrid,RealJ>
  :: operator() ( Field<Real,3,DGrid>& E,
                  Field<Real,3,DGrid>& B,
                  const Field<RealJ,3,DGrid>& J,
                  const mani::Grid<Real,DGrid>& grid,
                  const mpi::CartComm& comm,
                  Real dt,
                  Real preJ
                  ) const {
    // NOTE E,B,J are assumed to have been synced
    constexpr Real alpha = 0.5;

    const auto& mesh = E.mesh(); // TODO this implicitly requires J to have same grid as E
    apt::array<int,3> strd;
    for ( int i = 0; i < DGrid; ++i ) strd[i] = mesh.stride(i);
    Derivative<DGrid,Real> d{strd,grid};

    Field<Real,3,DGrid> B_old(mesh);
    std::swap(B,B_old);

    // update B
    for ( int C = 0; C < 3; ++C ) {
      Real abtt = dt*dt*alpha*(1.0-alpha);
      for ( int i = 0; i < B[0].data().size(); ++i ) {
        B[C][i] = B[C][i] + abtt * d.Lapl(B[C][i]);
      }
    }
    add_plus_curl<0>(B,-dt,E,d);
    add_plus_curl<1>(B,-dt,E,d);
    add_plus_curl<2>(B,-dt,E,d);

    add_plus_curl<0>(B,-alpha*dt*dt*preJ,J,d);
    add_plus_curl<1>(B,-alpha*dt*dt*preJ,J,d);
    add_plus_curl<2>(B,-alpha*dt*dt*preJ,J,d);

    copy_sync_guard_cells(B, comm);

    // update E partly before B_old is used and modified
    add_plus_curl<0>(E,(Real(1.0) - alpha) * dt,B_old,d);
    add_plus_curl<1>(E,(Real(1.0) - alpha) * dt,B_old,d);
    add_plus_curl<2>(E,(Real(1.0) - alpha) * dt,B_old,d);

    int num_iteration = 5;
    for ( int n = 0; n < num_iteration; ++n ) {
      std::swap(B,B_old);
      Real aa = alpha * alpha * dt * dt;
      for ( int C = 0; C < 3; ++C ) {
        for ( int i = 0; i < B[0].data().size(); ++i ) {
          B[C][i] = B_old[C][i] + aa * d.Lapl(B_old[C][i]);
        }
      }
      copy_sync_guard_cells(B, comm);
    }

    // update E, the rest of it
    add_plus_curl<0>(E,alpha*dt,B,d);
    add_plus_curl<1>(E,alpha*dt,B,d);
    add_plus_curl<2>(E,alpha*dt,B,d);
    for ( int C = 0; C < 3; ++C ) {
      for( int i = 0; i < E[0].data().size(); ++i )
        E[C][i] -= preJ * dt * J[C][i];
    }

    copy_sync_guard_cells(E, comm);
  }
}

#endif
