#ifndef _FIELD_CARTESIAN_SOLVER_IMPL_HPP_
#define _FIELD_CARTESIAN_SOLVER_IMPL_HPP_

#include "field/cartesian_solver/updater.hpp"
#include "field/sync.hpp"
#include "field/yee.hpp"

namespace field {
  // TODOL C++20 has templated lambda, following classes can be put this inside operator()
  template < int DGrid, typename T >
  struct Derivative {
  private:
    constexpr T diff_plus( const T& f, int stride ) const noexcept {
      return *(&f + stride) - f;
    }

    constexpr T diff_minus( const T& f, int stride ) const noexcept {
      return  f - *(&f - stride);
    }

    constexpr T diff2( const T& f, int stride ) const noexcept {
      return *(&f + stride) + *(&f - stride) - 2.0 * f;
    }

    template < int I, offset_t OFS_FIELD >
    constexpr T D( const T& f, int stride ) const noexcept {
      static_assert( I >= 0 && I < 3 );
      if constexpr ( I >= DGrid ) return 0.0;
      else if ( offset_t(OFS_FIELD) == INSITU ) { // f being of Etype
        return diff_plus(f,stride) / g[I].delta();
      } else { // f being of Btype
        return diff_minus(f,stride) / g[I].delta();
      }
    }

    template < int I >
    constexpr T DD( const T& f, int stride ) const noexcept {
      if constexpr ( I >= DGrid ) return 0.0;
      else return diff2( f, stride ) / (g[I].delta() * g[I].delta());
    }

  public:
    const mani::Grid<T,DGrid>& g; // grid

    template < int I >
    constexpr T curlcurl( Field<T,3,DGrid>& f, const apt::Index<DGrid>& idx ) const noexcept {
      const auto& m = f.mesh();
      int li = m.linearized_index_of_whole_mesh(idx);
      // used curlcurl == - laplacian for B field
      return - ( DD<0>(f[I][li],m.stride(0)) + DD<1>(f[I][li],m.stride(1)) + DD<2>(f[I][li],m.stride(2)) );
    }

    template < int I, offset_t Ftype >
    constexpr T curl( Field<T,3,DGrid>& f, const apt::Index<DGrid>& idx ) const noexcept { // idx is with respect to f
      // NOTE do not use f.offset. Use Ftype
      constexpr int J = (I+1)%3;
      constexpr int K = (I+2)%3;
      const auto& m = f.mesh();
      int li = m.linearized_index_of_whole_mesh(idx);
      // TODOL in c++20 there is templated lambda
      return D<J,Ftype>( f[K][li], m.stride(J) ) - D<K,Ftype>( f[J][li], m.stride(K) );
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
