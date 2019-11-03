#ifndef _FIELD_HAUGBOLLE_SOLVER_IMPL_HPP_
#define _FIELD_HAUGBOLLE_SOLVER_IMPL_HPP_

#include "field/yee.hpp"
#include "field/haugbolle_solver/updater.hpp"
#include "field/sync.hpp" // FIXME not in use right now

namespace field {
  template < int DGrid, typename T >
  struct VectorCalculus {
  private:
    static_assert( DGrid > 1 && DGrid < 4 );

  public:
    const apt::Grid<T,DGrid>& g; // grid
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
            return dir < DGrid ? g[dir].absc( idx[dir], 0.5 * ( (dir==f_comp) xor Ftype ) ) : 0;
          };
      return
        ( D[Ftype][K][J]( f[K][li], absc_fh(0,K), absc_fh(1,K), g, f.mesh().stride() )
          - D[Ftype][J][K]( f[J][li], absc_fh(0,J), absc_fh(1,J), g, f.mesh().stride() )
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
                  const apt::Grid<R,DGrid>& grid,
                  const mpi::CartComm& comm,
                  int timestep,
                  R dt
                  ) const {
    // NOTE E,B,J are assumed to have been synced
    const R alpha = this->implicit();
    const R beta = 1.0 - alpha;

    // running Begin and End
    auto rB = apt::range::far_begin(*this);
    auto rE = apt::range::far_end(*this);

    VectorCalculus<DGrid,R> vc{grid,_D,_hh};

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
            auto x = _hh[yee::Etype][comp](q1, q2, q3);
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
        B[C][li] -= beta * tmp[C][li]; // NOTE tmp has same mesh as B, so using li is OK here
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
    }

    // 5. finish updating B
    for ( int i = 0; i < DGrid; ++i ) ++rB[i]; // FIXME
    R prefactor = alpha*dt;
    for ( int C = 0; C < 3; ++C ) {
      for ( const auto& I : apt::Block(rB,rE) ) {
        B[C](I) -= prefactor * vc.template Curl<yee::Etype>(C,E,I);
      }
    }
  }

  template < typename R, int DGrid, typename RJ >
  void HaugbolleBdry<R,DGrid,RJ>
  :: operator() ( Field<R,3,DGrid>& E,
                  Field<R,3,DGrid>& B,
                  const Field<RJ,3,DGrid>& J,
                  const apt::Grid<R,DGrid>& grid,
                  const mpi::CartComm& comm,
                  int timestep,
                  R dt
                  ) const {
    // NOTE E,B,J are assumed to have been synced
    const R alpha = this->implicit();
    const R beta = 1.0 - alpha;

    // running Begin and End
    auto rB = apt::range::far_begin(*this);
    auto rE = apt::range::far_end(*this);
    VectorCalculus<DGrid,R> vc{grid};

    Field<R,3,DGrid> tmp{*this};
    Field<R,3,DGrid> tmp2{*this};
    for ( int C = 0; C < 3; ++C ) {
      for ( const auto& I : apt::Block(apt::range::far_begin(*this), apt::range::far_end(*this)) )
        tmp[C](I) = B[C](I); // store old B
    }
    tmp2.reset();

    auto setvc =
      [this](auto& vc,const auto& I, bool is_continuous_E = false) {
        int idx = _idx_von_neumann(I);
        for ( int i = 0; i < 3; ++i ) {
          vc.hh[yee::Etype][i] = _hh[yee::Etype][i][idx];
          vc.hh[yee::Btype][i] = _hh[yee::Btype][i][idx];
          for ( int s = 1; s < 3; ++s ) {
            vc.D[yee::Etype][(i+s)%3][i] = _D[yee::Etype][(i+s)%3][i][idx];
            vc.D[yee::Btype][(i+s)%3][i] = _D[yee::Btype][(i+s)%3][i][idx];
          }
        }
        if ( is_continuous_E ) {
          auto I1 = I;
          for ( int d = 0; d < DGrid; ++d ) {
            if ( 0 == _bdry[d] || !_continuous_transverse_E[d] ) continue;
            I1[d] = ( _bdry[d] == -1 ? apt::range::end(*this,d) : apt::range::begin(*this,d) );
            for ( int s = 1; s < 3; ++s )
              vc.D[yee::Etype][(d+s)%3][d] = _D[yee::Etype][(d+s)%3][d][_idx_von_neumann(I1)];
          }
        }
      };

    { // 1. find a convient field combination
      R prefactor = alpha*beta*dt;
      for ( int i = 0; i < DGrid; ++i ) {if ( _bdry[i] != -1 ) ++rB[i];}
      for ( const auto& I : apt::Block(rB,rE) ) {
        setvc(vc,I,true);
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
          [comp=C,this](R q1, R q2, R q3, const auto& I) -> R {
            int idx = _idx_von_neumann(I);
            auto x = _hh[yee::Etype][comp][idx](q1, q2, q3);
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
    { // 3. partial update B
      for ( int C = 0; C < 3; ++C ) {
        for ( const auto& I : apt::Block(rB,rE) ) {
          int li = B.mesh().linear_index(I);
          B[C][li] -= beta * tmp[C](I); //NOTE tmp has different mesh, cannot use li
          B[C][li] /= alpha;
        }
      }
    }
    { // 4. iteratively update E.
      R aatt = alpha*alpha*dt*dt;
      for ( int n = 0; n < _num_iter; ++n ) {
        for ( int i = 0; i < DGrid; ++i ) {if ( _bdry[i] != -1 ) ++rB[i];}
        if ( 0 == n ) {
          //  The first iteration is done differently to respect continuity of E and J
          for ( const auto& I : apt::Block(rB,rE) ) {
            setvc(vc,I,true);
            for ( int C = 0; C < 3; ++C )
              tmp[C](I) = vc.template Curl<yee::Etype>(C,E,I);
          }
          for ( const auto& I : apt::Block(rB,rE) ) {
            setvc(vc,I);
            for ( int C = 0; C < 3; ++C ) {
              tmp[C](I) += vc.template Curl<yee::Etype>(C,tmp2,I);
              E[C](I) += tmp2[C](I); // needed because below we have E -= ...
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
    }
    { // 5. finish updating B
      for ( int i = 0; i < DGrid; ++i ) ++rB[i]; // FIXME
      R prefactor = alpha*dt;
      for ( const auto& I : apt::Block(rB,rE) ) {
        setvc(vc,I);
        for ( int C = 0; C < 3; ++C )
          B[C](I) -= prefactor * vc.template Curl<yee::Etype>(C,E,I);
      }
    }
  }

}

#endif
