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

  template < int DGrid, typename Real, typename RealJ >
  struct RHS {
    const Mesh<DGrid>& mesh;
    const Field<RealJ,3,DGrid>& J;
    const Derivative<DGrid,Real>& d;

    template < int I_ >
    inline void ofB( Field<Real,3,DGrid>& B_new, Real dt, Real alpha, Real preJ, const Field<Real,3,DGrid>& E, const Field<Real,3,DGrid>& B ) const {
      constexpr int J_ = (I_+1)%3;
      constexpr int K_ = (I_+2)%3;
      for ( const auto& I : apt::Block(mesh.bulk_dims()) ) {
        int idx = mesh.linearized_index_of_whole_mesh(I);
        B_new[I_][idx] = B[I_][idx] + dt*dt*alpha*(1.0-alpha) * d.Lapl(B[I_][idx])
          - dt * ( d.template Curl<I_>(E[K_][idx - d.s[J_]], E[J_][idx - d.s[K_]] ) + alpha*dt*preJ* d.template Curl<I_>(J[K_][idx - d.s[J_]], J[J_][idx - d.s[K_]] ));
      }
    }

    template < int I_ >
    inline void ofE( Field<Real,3,DGrid>& E, Real dt, Real alpha, Real preJ, const Field<Real,3,DGrid>& B_old, const Field<Real,3,DGrid>& B_new ) const {
      constexpr int J_ = (I_+1)%3;
      constexpr int K_ = (I_+2)%3;
      for ( const auto& I : apt::Block(mesh.bulk_dims()) ) {
        int idx = mesh.linearized_index_of_whole_mesh(I);
        E[I_][idx] += dt * ( alpha * d.template Curl<I_>( B_new[K_][idx], B_new[J_][idx] ) + ( 1.0 - alpha ) * d.template Curl<I_>( B_old[K_][idx], B_old[J_][idx] ) - preJ * J[I_][idx] );
      }
    }

  };


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
    constexpr Real alpha = 0.5;

    const auto& mesh = E.mesh(); // TODO this implicitly requires J to have same grid as E
    apt::array<int,3> strd;
    for ( int i = 0; i < DGrid; ++i ) strd[i] = mesh.stride(i);
    Derivative<DGrid,Real> d{strd,grid};

    Field<Real,3,DGrid> tmp(mesh);

    RHS<DGrid,Real,RealJ> rhs{mesh,J,d};

    // update B
    rhs.template ofB<0>(tmp,dt,alpha,preJ,E,B);
    rhs.template ofB<1>(tmp,dt,alpha,preJ,E,B);
    rhs.template ofB<2>(tmp,dt,alpha,preJ,E,B);

    copy_sync_guard_cells(tmp, comm);

    // TODO solve a sparse matrix
    int num_iteration = 5;
    for ( int n = 0; n < num_iteration; ++n ) {
      
      // copy_sync_guard_cells(tmp, comm);
    }

    std::swap(B,tmp);

    // update E
    rhs.template ofE<0>(E,dt,alpha,preJ,tmp,B);
    rhs.template ofE<1>(E,dt,alpha,preJ,tmp,B);
    rhs.template ofE<2>(E,dt,alpha,preJ,tmp,B);

    copy_sync_guard_cells(E, comm);
  }
}

#endif
