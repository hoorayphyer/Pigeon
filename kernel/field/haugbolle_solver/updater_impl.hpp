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
    const apt::array<int,DGrid+1>& s; // stride
    apt::array<apt::array<apt::array<Deriv_t<T,DGrid>,3>,3>,2> D{}; // derivatives D[Ftype][Fcomp][coordinate]
    apt::array<apt::array<HH_t<T>,3>,2> hh;

    template < int I, offset_t Ftype >
    constexpr T curl( Field<T,3,DGrid>& f, const apt::Index<DGrid>& idx ) const noexcept { // idx is with respect to f
      // NOTE do not use f.offset. Use Ftype
      // NOTE Ftype is the offset for the ith-direction of f[i]. So transverse directions have !Ftype.
      constexpr int J = (I+1)%3;
      constexpr int K = (I+2)%3;
      int li = f.mesh().linear_index(idx);
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
      return
        ( D[Ftype][K][J]( f[K][li], absc_fh(0,K), absc_fh(1,K), g, s )
          - D[Ftype][J][K]( f[J][li], absc_fh(0,J), absc_fh(1,J), g, s )
          ) / hh[!Ftype][I]( absc_hh(0,I), absc_hh(1,I), absc_hh(2,I) );
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

  template < typename R, int DGrid, typename RJ >
  void Haugbolle<R,DGrid,RJ>
  :: operator() ( Field<R,3,DGrid>& E,
                  Field<R,3,DGrid>& B,
                  const Field<RJ,3,DGrid>& J,
                  const mani::Grid<R,DGrid>& grid,
                  const mpi::CartComm& comm,
                  int timestep,
                  R dt
                  ) const {
    // NOTE E,B,J are assumed to have been synced
    constexpr R alpha = 0.501;
    constexpr R beta = 1.0 - alpha;

    // running Begin and End
    auto rB = apt::range::far_begin(*this);
    auto rE = apt::range::far_end(*this);

    VectorCalculus<DGrid,R> vc{grid,E.mesh().stride(),_D,_hh};

    Field<R,3,DGrid> tmp(B); // tmp holds old B values

    { // 1. find a convient field combination
      R prefactor = alpha*beta*dt;
      for ( int i = 0; i < DGrid; ++i ) ++rB[i];
      for ( int C = 0; C < 3; ++C ) {
        for ( const auto& I : apt::Block(rB,rE) ) {
          B[C](I) -= prefactor * vc.template Curl<yee::Etype>(C,E,I);
        }
      }
    }
    { // 2. partial update E.
      for ( int i = 0; i < DGrid; ++i ) --rE[i];
      for ( int C = 0; C < 3; ++C ) {
        for ( const auto& I : apt::Block(rB,rE) ) {
          E[C](I) += dt * vc.template Curl<yee::Btype>(C,B,I);
        }
      }

      // subtract J from E. Divide J by hh
      R prefactor = dt*_preJ_factor;
      for ( int C = 0; C < 3; ++C ) {
        auto hhinv =
          [comp=C,this](R q1, R q2, R q3) -> R {
            auto x = this->_hh[yee::Etype][comp](q1, q2, q3);
            if ( std::abs(x) < 1e-8 ) return 0;
            else return 1 / x;
          };
        const auto ofs = J[C].offset();
        R q[3] = {0, 0, 0};
        for ( const auto& I : apt::Block(rB,rE) ) {
          for ( int i = 0; i < DGrid; ++i ) q[i] = grid[i].absc(I[i], 0.5 * ofs[i]);
          E[C](I) -= ( prefactor * J[C](I) * hhinv(q[0],q[1],q[2]) );
        }
      }
    }
    // 3. partial update B.
    for ( int C = 0; C < 3; ++C ) {
      for ( const auto& I : apt::Block(rB, rE) ) {
        int li = B.mesh().linear_index(I);
        B[C][li] -= beta * tmp[C][li];
        B[C][li] /= alpha;
      }
    }
    { // 4. iteratively update F1. F1_new = F1_old - (alpha*dt)^2 * curl ( curl F1_old ). Iterate as follows. Let X Y be two field objects, X stores F1_old to start with
      //   - curl X and stores into Y.
      //   - X[i] -= (alpha * dt)^2 * (curl Y) [i]
      // now X has the updated F1
      R aatt = alpha*alpha*dt*dt;
      for ( int n = 0; n < _num_iter; ++n ) {

        for ( int i = 0; i < DGrid; ++i ) ++rB[i];
        for ( int C = 0; C < 3; ++C ) {
          for ( const auto& I : apt::Block(rB,rE) ) {
            tmp[C](I) = vc.template Curl<yee::Etype>(C,E,I);
          }
        }
        for ( int i = 0; i < DGrid; ++i ) --rE[i];
        for ( int C = 0; C < 3; ++C ) {
          for ( const auto& I : apt::Block(rB,rE) ) {
            E[C](I) -= aatt * vc.template Curl<yee::Btype>(C,tmp,I);
          }
        }

      }
      copy_sync_guard_cells(E, comm); // NOTE this needs to be done before updating B
    }

    // 5. finish updating B
    R prefactor = alpha*dt;
    for ( int C = 0; C < 3; ++C ) {
      for ( const auto& I : apt::Block(rB,rE) ) {
        B[C](I) -= prefactor * vc.template Curl<yee::Etype>(C,E,I);
      }
    }
    copy_sync_guard_cells(B, comm);
  }

  template < typename R, int DGrid, typename RJ >
  void HaugbolleBdry<R,DGrid,RJ>
  :: operator() ( Field<R,3,DGrid>& E,
                  Field<R,3,DGrid>& B,
                  const Field<RJ,3,DGrid>& J,
                  const mani::Grid<R,DGrid>& grid,
                  const mpi::CartComm& comm,
                  int timestep,
                  R dt
                  ) const {
    // NOTE E,B,J are assumed to have been synced
    constexpr R alpha = 0.501;
    constexpr R beta = 1.0 - alpha;

    apt::array<int,DGrid+1> storage_stride;
    {
      storage_stride[0] = 1;
      for ( int i = 0; i < DGrid; ++i )
        storage_stride[i+1] = storage_stride[i] * ( (0 == _bdry[i])  ? 1 : (apt::range::size(*this, i)) );
      for ( int i = 0; i < DGrid; ++i )
        if ( 0 == _bdry[i] ) storage_stride[i] = 0; // mask out nonbdry dimensions
    }

    auto storage_idx =
      [&s=storage_stride, &r=*this]( const auto& I ) {
        // NOTE this function should encapsulate all details of resolving index
        int idx = 0;
        for ( int i = 0; i < DGrid; ++i )
          idx += std::max( 0, std::min( I[i] - apt::range::begin(r,i), apt::range::size(r,i) - 1 ) ) * s[i];
        return idx;
      };

    auto index_of_continuous_diff_in_direction =
      [&s=storage_stride, this]( int dir, const auto& I ) {
        int idx = 0;
        for ( int i = 0; i < DGrid; ++i ) {
          if ( i != dir || this->_bdry[i] == 0 )
            idx += std::max( 0, std::min( I[i] - apt::range::begin(*this,i), apt::range::size(*this,i) - 1 ) ) * s[i];
          else {
            idx += ( this->_bdry[i] == -1 ? apt::range::size(*this,i) - 1 : 0 ) * s[i];
          }
        }
        return idx;
      };

    auto make_continuous =
      [&, this](auto& vc, const auto& I) {
        for ( int i = 0; i < DGrid; ++i ) {
          if ( 0 == this->_bdry[i] || !this->_continuous_transverse_E[i] ) continue;
          for ( int s = 1; s < 3; ++s )
            vc.D[yee::Etype][(i+s)%3][i] = this->_D[yee::Etype][(i+s)%3][i][index_of_continuous_diff_in_direction(i,I)];
        }
      };

    // running Begin and End
    auto rB = apt::range::far_begin(*this);
    auto rE = apt::range::far_end(*this);
    VectorCalculus<DGrid,R> vc{grid,E.mesh().stride()};

    Field<R,3,DGrid> tmp{*this};
    Field<R,3,DGrid> tmp2{*this};
    for ( int C = 0; C < 3; ++C ) {
      for ( const auto& I : apt::Block(apt::range::far_begin(*this), apt::range::far_end(*this)) )
        tmp[C](I) = B[C](I); // store old B
    }
    tmp2.reset();

    auto setvc =
      [this, &storage_idx](auto& vc,const auto& I) {
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
      R prefactor = alpha*beta*dt;
      for ( int i = 0; i < DGrid; ++i ) {if ( _bdry[i] != -1 ) ++rB[i];}
      for ( const auto& I : apt::Block(rB,rE) ) {
        setvc(vc,I);
        make_continuous(vc,I);
        for ( int C = 0; C < 3; ++C )
          B[C](I) -= prefactor * vc.template Curl<yee::Etype>(C,E,I);
      }
    }
    { // 2. Partial Update on E (different from Haugbolle bulk treatment)
      // store curlB separately because it may be discontinuous across interface
      for ( int i = 0; i < DGrid; ++i ) {if ( _bdry[i] != 1 ) --rE[i];}
      for ( const auto& I : apt::Block(rB,rE) ) {
        setvc(vc,I);
        for ( int C = 0; C < 3; ++C )
          tmp2[C](I) = dt * vc.template Curl<yee::Btype>(C,B,I);
      }

      // subtract J from E. Divide J by hh.
      R prefactor = dt*_preJ_factor;
      for ( int C = 0; C < 3; ++C ) {
        auto hhinv =
          [comp=C,this,&storage_idx](R q1, R q2, R q3, const auto& I) -> R {
            auto x = this->_hh[yee::Etype][comp][storage_idx(I)](q1, q2, q3);
            if ( std::abs(x) < 1e-8 ) return 0;
            else return 1 / x;
          };
        const auto ofs = J[C].offset();
        R q[3] = {0, 0, 0};
        for ( const auto& I : apt::Block(rB,rE) ) {
          for ( int i = 0; i < DGrid; ++i ) q[i] = grid[i].absc(I[i], 0.5 * ofs[i]);
          E[C](I) -= ( prefactor * J[C](I) * hhinv(q[0],q[1],q[2],I) );
        }
      }
    }
    // 3. partial update B
    for ( int C = 0; C < 3; ++C ) {
      for ( const auto& I : apt::Block(rB,rE) ) {
        int li = B.mesh().linear_index(I);
        B[C][li] -= beta * tmp[C][li];
        B[C][li] /= alpha;
      }
    }
    { // 4. iteratively update E.
      R aatt = alpha*alpha*dt*dt;
      for ( int n = 0; n < _num_iter; ++n ) {
        for ( int i = 0; i < DGrid; ++i ) {if ( _bdry[i] != -1 ) ++rB[i];}
        if ( 0 == n ) {
          //  The first iteration is done differently to respect continuity of E and J
          for ( const auto& I : apt::Block(rB,rE) ) {
            setvc(vc,I);
            make_continuous(vc,I);
            for ( int C = 0; C < 3; ++C )
              tmp[C](I) = vc.template Curl<yee::Etype>(C,E,I);
          }
          for ( const auto& I : apt::Block(rB,rE) ) {
            setvc(vc,I);
            for ( int C = 0; C < 3; ++C ) {
              tmp[C](I) += vc.template Curl<yee::Etype>(C,tmp2,I);
              E[C](I) -= tmp2[C](I);
            }
          }
        } else {
          for ( const auto& I : apt::Block(rB,rE) ) {
            setvc(vc,I);
            for ( int C = 0; C < 3; ++C )
              tmp[C](I) = vc.template Curl<yee::Etype>(C,E,I);
          }
        }

        for ( int i = 0; i < DGrid; ++i ) {if ( _bdry[i] != 1 ) --rE[i];}
        for ( const auto& I : apt::Block(rB,rE) ) {
          setvc(vc,I);
          for ( int C = 0; C < 3; ++C )
            E[C](I) -= aatt * vc.template Curl<yee::Btype>(C,tmp,I);
        }
      }
      copy_sync_guard_cells(E, comm, E.mesh().range(), *this);
    }
    // 5. finish updating B
    R prefactor = alpha*dt;
    for ( const auto& I : apt::Block(rB,rE) ) {
      setvc(vc,I);
      for ( int C = 0; C < 3; ++C )
        B[C](I) -= prefactor * vc.template Curl<yee::Etype>(C,E,I);
    }

    copy_sync_guard_cells(B, comm, B.mesh().range(), *this);
  }
}

#endif
