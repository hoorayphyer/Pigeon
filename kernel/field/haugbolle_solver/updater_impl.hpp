#ifndef _FIELD_HAUGBOLLE_SOLVER_IMPL_HPP_
#define _FIELD_HAUGBOLLE_SOLVER_IMPL_HPP_

#include "field/yee.hpp"
#include "field/haugbolle_solver/updater.hpp"
#include "field/sync.hpp"

namespace field {
  template < int DGrid, typename T >
  struct VectorCalculus {
  private:
    static_assert( DGrid > 1 && DGrid < 4 );

  public:
    const mani::Grid<T,DGrid>& g; // grid
    const apt::Index<DGrid> s; // stride
    apt::array<apt::array<apt::array<Deriv_t<T,DGrid>,3>,3>,2> D{}; // derivatives D[Ftype][Fcomp][coordinate]
    apt::array<apt::array<HH_t<T>,3>,2> hh;

    template < int I, offset_t Ftype >
    constexpr T curl( Field<T,3,DGrid>& f, const apt::Index<DGrid>& idx ) const noexcept { // idx is with respect to f
      // NOTE do not use f.offset. Use Ftype
      // NOTE Ftype is the offset for the ith-direction of f[i]. So transverse directions have !Ftype.
      constexpr int J = (I+1)%3;
      constexpr int K = (I+2)%3;
      int li = f.mesh().linearized_index_of_whole_mesh(idx);
      // TODOL in c++20 there is templated lambda
      auto absc_fh
        = [&idx,&g=this->g]( int dir, int f_comp ) noexcept {
            // needs f's offset, which is Ftype when dir == f_comp
            return dir < DGrid ? g[dir].absc( idx[dir], 0.5 * ( (dir==f_comp) xor !Ftype ) ) : 0;
          };
      auto absc_hh
        = [&idx,&g=this->g]( int dir, int f_comp ) noexcept {
            // needs f's anti-offset, which is !Ftype when dim == f_comp
            return dir < DGrid ? g[dir].absc( idx[dir], 0.5 * ( (dir==f_comp) xor Ftype ) ) : 0.0;
          };
      // FGH
      return
        ( D[Ftype][K][J]( f[K][li], absc_fh(0,K), absc_fh(1,K), g, s )
          - D[Ftype][J][K]( f[J][li], absc_fh(0,J), absc_fh(1,J), g, s )
          );// / hh[!Ftype][I]( absc_hh(0,I), absc_hh(1,I), absc_hh(2,I) );
    }

    template < offset_t Ftype >
    constexpr T Curl(int i, Field<T,3,DGrid>& f, const apt::Index<DGrid>& idx ) const noexcept {
      switch (i) {
      case 0 : return curl<0,Ftype>(f,idx);
      case 1 : return curl<1,Ftype>(f,idx);
      case 2 : return curl<2,Ftype>(f,idx);
      }
    }

  };

  template < typename Real, int DGrid, typename RealJ >
  void Haugbolle<Real,DGrid,RealJ>
  :: operator() ( Field<Real,3,DGrid>& E,
                  Field<Real,3,DGrid>& B,
                  const Field<RealJ,3,DGrid>& J,
                  const mani::Grid<Real,DGrid>& grid,
                  const mpi::CartComm& comm,
                  int timestep,
                  Real dt
                  ) const {
    // NOTE E,B,J are assumed to have been synced
    constexpr Real alpha = 0.5;
    constexpr Real beta = 1.0 - alpha;

    const auto block = apt::Block(this->ext());
    VectorCalculus<DGrid,Real> vc{grid,E.mesh().strides(),_D,_hh};

    // adjust ranges for first_curl
    const auto [Ib_1st_curl, block_1st_curl ] =
      []( auto Ib, auto ext ) {
        for ( int i = 0; i < DGrid; ++i ) {
          Ib[i] -= 1;
          ext[i] += 2;
        }
        return std::make_pair( Ib, apt::Block(ext));
      } (this->Ib(), this->ext());

    Field<Real,3,DGrid> tmp(B); // tmp holds old B values

    { // 1. find a convient field combination
      Real prefactor = alpha*beta*dt;
      for ( int C = 0; C < 3; ++C ) {
        for ( auto I : block_1st_curl ) {
          I += Ib_1st_curl;
          B[C](I) -= prefactor * vc.template Curl<yee::Etype>(C,E,I);
        }
      }
    }
    { // 2. partial update E. Note the ranges are only the bulk
      for ( int C = 0; C < 3; ++C ) {
        for ( auto I : block ) {
          I += this->Ib();
          E[C](I) += dt * vc.template Curl<yee::Btype>(C,B,I);
        }
      }

      // subtract J from E
      Real prefactor = dt*_preJ_factor;
      for ( int C = 0; C < 3; ++C ) {
        for ( int i = 0; i < E[0].data().size(); ++i ) {
          E[C][i] -= prefactor * J[C][i];
        }
      }
      copy_sync_guard_cells(E, comm);
    }
    // 3. partial update B
    for ( int C = 0; C < 3; ++C ) {
      for ( int i = 0; i < E[0].data().size(); ++i ) {
        B[C][i] -= beta * tmp[C][i];
        B[C][i] /= alpha;
      }
    }
    { // 4. iteratively update F1. F1_new = F1_old - (alpha*dt)^2 * curl ( curl F1_old ). Iterate as follows. Let X Y be two field objects, X stores F1_old to start with
      //   - curl X and stores into Y. NOTE curl X should be done from 1 cell into the guard
      //   - X[i] -= (alpha * dt)^2 * (curl Y) [i]
      //   - sync X
      // now X has the updated F1
      Real aatt = alpha*alpha*dt*dt;
      constexpr int num_iteration = 4;
      for ( int n = 0; n < num_iteration; ++n ) {
        for ( int C = 0; C < 3; ++C ) {
          for ( auto I : block_1st_curl ) {
            I += Ib_1st_curl;
            tmp[C](I) = vc.template Curl<yee::Etype>(C,E,I);
          }
        }
        for ( int C = 0; C < 3; ++C ) {
          for ( auto I : block ) {// NOTE range is only in bulk
            I += this->Ib();
            E[C](I) -= aatt * vc.template Curl<yee::Btype>(C,tmp,I);
          }
        }
      }
      copy_sync_guard_cells(E, comm);
    }
    // 5. finish updating B
    Real prefactor = alpha*dt;
    for ( int C = 0; C < 3; ++C ) {
      for ( auto I : block ) {
        I += this->Ib();
        B[C](I) -= prefactor * vc.template Curl<yee::Etype>(C,E,I);
      }
    }
    copy_sync_guard_cells(B, comm);
  }

  template < typename Real, int DGrid, typename RealJ >
  void HaugbolleBdry<Real,DGrid,RealJ>
  :: operator() ( Field<Real,3,DGrid>& E,
                  Field<Real,3,DGrid>& B,
                  const Field<RealJ,3,DGrid>& J,
                  const mani::Grid<Real,DGrid>& grid,
                  const mpi::CartComm& comm,
                  int timestep,
                  Real dt
                  ) const {
    // NOTE E,B,J are assumed to have been synced
    constexpr Real alpha = 0.5;
    constexpr Real beta = 1.0 - alpha;

    static apt::array<int,DGrid+1> storage_stride;
    {
      storage_stride[0] = 1;
      for ( int i = 0; i < DGrid; ++i )
        storage_stride[i+1] = storage_stride[i] * ( (0 == _bdry[i])  ? 1 : this->ext()[i] );
      for ( int i = 0; i < DGrid; ++i )
        if ( 0 == _bdry[i] ) storage_stride[i] = 0; // mask out nonbdry dimensions
    }

    auto storage_idx =
      [&s=storage_stride, &Ib=this->Ib(), ext=this->ext()]( const auto& I ) {
        int idx = 0;
        for ( int i = 0; i < DGrid; ++i )
          idx += std::max( 0, std::min( I[i] - Ib[i], ext[i] - 1 ) ) * s[i];
        return idx;
      };

    const auto block = apt::Block(this->ext());
    VectorCalculus<DGrid,Real> vc{grid,E.mesh().strides()};

    // adjust ranges for first_curl
    const auto [Ib_1st_curl, block_1st_curl ] =
      [&bdry=_bdry]( auto Ib, auto ext ) {
        for ( int i = 0; i < DGrid; ++i ) {
          if (bdry[i] != -1) {
            Ib[i] -= 1;
            ext[i] += 1;
          }
          if (bdry[i] != 1) ext[i] += 1;
        }
        return std::make_pair( Ib, apt::Block(ext));
      } (this->Ib(), this->ext());

    Field<Real,3,DGrid> tmp(B); // tmp holds old B values

    auto setvc =
      [this, storage_idx](auto& vc,const auto& I) {
        int idx = storage_idx(I);
        for ( int i = 0; i < 3; ++i ) {
          vc.hh[yee::Etype][i] = _hh[yee::Etype][i][idx];
          vc.hh[yee::Btype][i] = _hh[yee::Btype][i][idx];
          for ( int s = 1; s < 3; ++s ) {
            vc.D[yee::Etype][(i+s)%3][i] = _D[yee::Etype][(i+s)%3][i][idx];
            vc.D[yee::Btype][(i+s)%3][i] = _D[yee::Btype][(i+s)%3][i][idx];
          }
        }
      };

    { // 1. find a convient field combination
      Real prefactor = alpha*beta*dt;
      for ( auto I : block_1st_curl ) {
        I += Ib_1st_curl;
        setvc(vc,I);
        // for the first curl of E, enforce continuity
        for ( int i = 0; i < DGrid; ++i ) {
          if ( 0 == _bdry[i] || !_continuous_transverse_E[i] ) continue;
          if ( -1 == _bdry[i] ) {
            for ( int s = 1; s < 3; ++s )
              // TODO FIXME not back but the back in one dimension
              vc.D[yee::Etype][(i+s)%3][i] = _D[yee::Etype][(i+s)%3][i].back();
          } else {
            for ( int s = 1; s < 3; ++s )
              vc.D[yee::Etype][(i+s)%3][i] = _D[yee::Etype][(i+s)%3][i].front();
          }
        }
        for ( int C = 0; C < 3; ++C )
          B[C](I) -= prefactor * vc.template Curl<yee::Etype>(C,E,I);
      }
    }
    { // 2. partial update E. Note the ranges are only the bulk
      for ( auto I : block ) {
        I += this->Ib();
        setvc(vc,I);
        for ( int C = 0; C < 3; ++C )
          E[C](I) += dt * vc.template Curl<yee::Btype>(C,B,I);
      }

      // subtract J from E
      Real prefactor = dt*_preJ_factor;
      for ( int C = 0; C < 3; ++C ) {
        for ( int i = 0; i < E[0].data().size(); ++i ) {
          E[C][i] -= prefactor * J[C][i];
        }
      }
      copy_sync_guard_cells(E, comm);
    }
    // 3. partial update B
    for ( int C = 0; C < 3; ++C ) {
      for ( int i = 0; i < E[0].data().size(); ++i ) {
        B[C][i] -= beta * tmp[C][i];
        B[C][i] /= alpha;
      }
    }
    { // 4. iteratively update F1. F1_new = F1_old - (alpha*dt)^2 * curl ( curl F1_old ). Iterate as follows. Let X Y be two field objects, X stores F1_old to start with
      //   - curl X and stores into Y. NOTE curl X should be done from 1 cell into the guard
      //   - X[i] -= (alpha * dt)^2 * (curl Y) [i]
      //   - sync X
      // now X has the updated F1
      Real aatt = alpha*alpha*dt*dt;
      constexpr int num_iteration = 4;
      for ( int n = 0; n < num_iteration; ++n ) {
        for ( auto I : block_1st_curl ) {
          I += Ib_1st_curl;
          setvc(vc,I);
          for ( int C = 0; C < 3; ++C )
            tmp[C](I) = vc.template Curl<yee::Etype>(C,E,I);
        }
        // NOTE range is only in bulk
        for ( auto I : block ) {
          I += this->Ib();
          setvc(vc,I);
          for ( int C = 0; C < 3; ++C )
            E[C](I) -= aatt * vc.template Curl<yee::Btype>(C,tmp,I);
        }
      }
      copy_sync_guard_cells(E, comm);
    }
    // 5. finish updating B
    Real prefactor = alpha*dt;
    for ( auto I : block ) {
      I += this->Ib();
      setvc(vc,I);
      for ( int C = 0; C < 3; ++C )
        B[C](I) -= prefactor * vc.template Curl<yee::Etype>(C,E,I);
    }

    copy_sync_guard_cells(B, comm);
  }
}

#endif
