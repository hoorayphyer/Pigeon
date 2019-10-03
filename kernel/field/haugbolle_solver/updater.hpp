#ifndef _FIELD_HAUGBOLLE_UPDATER_HPP_
#define _FIELD_HAUGBOLLE_UPDATER_HPP_

#include "field/action.hpp"

namespace field {
  template < typename Real, int DGrid >
  using Deriv_t = Real(*)(const Real& , Real, Real, const mani::Grid<Real,DGrid>&, const apt::Index<DGrid>&);

  template < typename T >
  using HH_t = T(*)(T,T,T);

  template < typename Real, int DGrid, typename RealJ >
  struct Haugbolle: public Action<Real,DGrid,RealJ> {
  private:
    Real _preJ_factor {};
    apt::array<apt::array<apt::array<Deriv_t<Real,DGrid>,3>,3>,2> _D{}; // derivatives D[Ftype][Fcomp][coordinate]
    apt::array<apt::array<HH_t<Real>,3>,2> _hh {}; // hh[Ftype][Fcomp][coordinate], where F is the field being curled. hh_ij = hh_ji

  public:
    constexpr Haugbolle& set_preJ( Real preJ ) noexcept {
      _preJ_factor = preJ;
      return *this;
    }
    constexpr Haugbolle& set_D( offset_t Ftype, int field_comp, int drv, Deriv_t<Real,DGrid> f ) noexcept {
      _D[Ftype][field_comp][drv] = f;
      return *this;
    }

    constexpr Haugbolle& set_hh( offset_t Ftype, int k, HH_t<Real> f ) noexcept {
      _hh[Ftype][k] = f;
      return *this;
    }
    constexpr Haugbolle& set_hh( offset_t Ftype, int field_comp, int drv, HH_t<Real> f ) noexcept {
      return set_hh(!Ftype, 3 - field_comp - drv, f);
    }
    constexpr const auto& hh() const { return _hh; }
    constexpr const auto& D() const { return _D; }
    constexpr const auto& preJ() const { return _preJ_factor; }

    virtual Haugbolle* Clone() const { return new Haugbolle(*this); }

    virtual void operator() ( Field<Real,3,DGrid>& E,
                              Field<Real,3,DGrid>& B,
                              const Field<RealJ,3,DGrid>& Jmesh,
                              const mani::Grid<Real,DGrid>& grid,
                              const mpi::CartComm& cart,
                              int timestep,
                              Real dt
                              ) const override;
  };

  template < typename Real, int DGrid, typename RealJ >
  struct HaugbolleBdry: public Action<Real,DGrid,RealJ> {
  private:
    Real _preJ_factor {};
    apt::array<apt::array<apt::array<std::vector<Deriv_t<Real,DGrid>>,3>,3>,2> _D{}; // derivatives D[Ftype][Fcomp][coordinate]
    apt::array<apt::array<std::vector<HH_t<Real>>,3>,2> _hh {}; // hh[Ftype][Fcomp]

    apt::array<int,DGrid> _bdry{}; // -1 means lower, 0 means not a boundary, 1 means upper
    apt::array<bool,DGrid> _continuous_transverse_E{};

  public:
    HaugbolleBdry() = default;
    HaugbolleBdry(const HaugbolleBdry& ) = default;

    void init_from( const Haugbolle<Real,DGrid,RealJ>& haug, int size ) {
      _preJ_factor = haug.preJ();
      for ( int Ftype = 0; Ftype < 2; ++Ftype ) {
        for ( int i = 0; i < 3; ++i ) { // i is coordinate
          _hh[Ftype][i] = {size, haug.hh()[Ftype][i]};
          for ( int j = 0; j < 3; ++j ) {
            if ( i == j ) continue;
            _D[Ftype][j][i] = {size, haug.D()[Ftype][j][i]};
          }
        }
      }
    }

    constexpr HaugbolleBdry& set_preJ( Real preJ ) noexcept {
      _preJ_factor = preJ;
      return *this;
    }
    constexpr HaugbolleBdry& set_D( offset_t Ftype, int field_comp, int drv, Deriv_t<Real,DGrid> f, int idx ) noexcept {
      if ( idx >= _D[Ftype][field_comp][drv].size() )
        _D[Ftype][field_comp][drv].resize(idx+1);

      _D[Ftype][field_comp][drv][idx] = f;
      return *this;
    }

    constexpr HaugbolleBdry& set_hh( offset_t Ftype, int k, HH_t<Real> f, int idx ) noexcept {
      if ( idx >= _hh[Ftype][k].size() )
        _hh[Ftype][k].resize(idx+1);
      _hh[Ftype][k][idx] = f;
      return *this;
    }
    constexpr HaugbolleBdry& set_hh( offset_t Ftype, int field_comp, int drv, HH_t<Real> f, int idx ) noexcept {
      return set_hh(!Ftype, 3 - field_comp - drv, f, idx);
    }
    constexpr HaugbolleBdry& set_boundary( apt::array<int,DGrid> val ) noexcept {
      _bdry = val;
      return *this;
    }
    constexpr HaugbolleBdry& set_boundary( int b1, int b2=0, int b3=0 ) noexcept {
      int b [3] = {b1, b2, b3};
      for ( int i = 0; i < DGrid; ++i ) _bdry[i] = b[i];
      return *this;
    }
    constexpr HaugbolleBdry& enforce_continuous_transverse_E( apt::array<bool,DGrid> val ) noexcept {
      _continuous_transverse_E = val;
      return *this;
    }
    constexpr HaugbolleBdry& enforce_continuous_transverse_E( bool b1, bool b2=false, bool b3=false ) noexcept {
      bool b [3] = {b1, b2, b3};
      for ( int i = 0; i < DGrid; ++i ) _continuous_transverse_E[i] = b[i];
      return *this;
    }

    constexpr const auto& hh() const { return _hh; }
    constexpr const auto& D() const { return _D; }
    constexpr const auto& preJ() const { return _preJ_factor; }
    constexpr const auto& boundary() const { return _bdry; }
    constexpr const auto& cont_trans_E() const { return _continuous_transverse_E; }

    virtual HaugbolleBdry* Clone() const { return new HaugbolleBdry(*this); }

    virtual void operator() ( Field<Real,3,DGrid>& E,
                              Field<Real,3,DGrid>& B,
                              const Field<RealJ,3,DGrid>& Jmesh,
                              const mani::Grid<Real,DGrid>& grid,
                              const mpi::CartComm& cart,
                              int timestep,
                              Real dt
                              ) const override;
  };

}

#endif
