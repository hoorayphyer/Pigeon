#ifndef _FIELD_HAUGBOLLE_SOLVER_IMPL_HPP_
#define _FIELD_HAUGBOLLE_SOLVER_IMPL_HPP_

#include "field/haugbolle_solver/derivative_cartesian.hpp"
#include "field/haugbolle_solver/updater.hpp"
#include "field/sync.hpp"

namespace field {
  template < typename Real, int DGrid, typename RealJ >
  void HaugbolleUpdater<Real,DGrid,RealJ>
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
    constexpr Real beta = 1.0 - alpha;

    const auto& mesh = E.mesh(); // TODO this implicitly requires J to have same grid as E
    const auto block = apt::Block(mesh.bulk_dims());
    Derivative<DGrid,Real> d{grid};

    Field<Real,3,DGrid> tmp(E); // tmp holds old E values

    { // 1. find a convient field combination
      // subtract J from E
      for ( int C = 0; C < 3; ++C ) {
        Real prefactor = alpha*dt*preJ;
        for ( int i = 0; i < E[0].data().size(); ++i ) {
          E[C][i] -= prefactor * J[C][i];
        }
      }
      // add curl B to E
      apt::Index<DGrid> Ib;
      apt::Index<DGrid> ext = mesh.bulk_dims();
      for ( auto& x : Ib ) x = -1; // NOTE the ranges are one cell into the guard
      for ( auto& x : ext ) x += 2;
      Real prefactor = alpha*beta*dt;
      // TODO c++20 has templated lambda
      for ( auto I : apt::Block(ext) ) {
        I += Ib;
        E[0](I) += prefactor * d.template curl<0,yee::Btype>(B,I);
      }
      for ( auto I : apt::Block(ext) ) {
        I += Ib;
        E[1](I) += prefactor * d.template curl<1,yee::Btype>(B,I);
      }
      for ( auto I : apt::Block(ext) ) {
        I += Ib;
        E[2](I) += prefactor * d.template curl<2,yee::Btype>(B,I);
      }
    }
    // 2. partial update B. Note the ranges are only the bulk
    for ( const auto& I : block ) {
      B[0](I) -= dt * d.template curl<0,yee::Etype>(E,I);
    }
    for ( const auto& I : block ) {
      B[1](I) -= dt * d.template curl<1,yee::Etype>(E,I);
    }
    for ( const auto& I : block ) {
      B[2](I) -= dt * d.template curl<2,yee::Etype>(E,I);
    }
    copy_sync_guard_cells(B, comm);
    // 3. partial update E.
    for ( int C = 0; C < 3; ++C ) {
      for ( int i = 0; i < E[0].data().size(); ++i ) {
        E[C][i] -= beta * tmp[C][i];
        E[C][i] /= alpha;
      }
    }
    { // 4. iteratively update B
      constexpr int num_iteration = 5;
      for ( int n = 0; n < num_iteration; ++n ) {
        std::swap(B,tmp);
        Real aatt = alpha * alpha * dt * dt;
        // NOTE range is only in bulk
        for ( const auto& I : block ) {
          B[0](I) = tmp[0](I) - aatt * d.template curlcurl<0>(tmp,I);
        }
        for ( const auto& I : block ) {
          B[1](I) = tmp[1](I) - aatt * d.template curlcurl<1>(tmp,I);
        }
        for ( const auto& I : block ) {
          B[2](I) = tmp[2](I) - aatt * d.template curlcurl<2>(tmp,I);
        }
        copy_sync_guard_cells(B, comm);
      }
    }

    // 5. finish updating E
    Real prefactor = alpha*dt;
    for ( const auto& I : block ) {
      E[0](I) += prefactor * d.template curl<0,yee::Btype>(B,I);
    }
    for ( const auto& I : block ) {
      E[1](I) += prefactor * d.template curl<1,yee::Btype>(B,I);
    }
    for ( const auto& I : block ) {
      E[2](I) += prefactor * d.template curl<2,yee::Btype>(B,I);
    }

    copy_sync_guard_cells(E, comm);
  }
}

#endif
