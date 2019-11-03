#ifndef _FIELD_HAUGBOLLE_UPDATER_HPP_
#define _FIELD_HAUGBOLLE_UPDATER_HPP_

#include "field/action.hpp"

namespace field {
  template < class H, typename R >
  struct HaugbolleParams {
  private:
    R _fourpi {};
    int _num_iter = 4;
    R _implicit = 0.501; // NOTE at > 0.5 to be stable

  public:
    constexpr H& set_fourpi( R x ) noexcept {
      _fourpi = x;
      return static_cast<H&>(*this);
    }
    constexpr H& set_number_iteration( int n ) noexcept {
      _num_iter = n;
      return static_cast<H&>(*this);
    }
    constexpr H& set_implicit( R a ) noexcept {
      _implicit = a;
      return static_cast<H&>(*this);
    }

    constexpr const auto& fourpi() const { return _fourpi; }
    constexpr const auto& num_iter() const { return _num_iter; }
    constexpr const auto& implicit() const { return _implicit; }
  };

  template < typename R, int DGrid >
  using Deriv_t = R(*)(const R& , R, R, const apt::Grid<R,DGrid>&, const apt::Index<DGrid+1>&);

  template < typename T >
  using HH_t = T(*)(T,T,T);

  template < typename R, int DGrid, typename RJ >
  struct Haugbolle: public HaugbolleParams<Haugbolle<R,DGrid,RJ>,R>,
                    public Action<R,DGrid,RJ> {
  private:
    R _preJ_factor {};
    apt::array<apt::array<apt::array<Deriv_t<R,DGrid>,3>,3>,2> _D{}; // derivatives D[Ftype][Fcomp][coordinate]
    apt::array<apt::array<HH_t<R>,3>,2> _hh {};
    int _num_iter = 4;

  public:
    constexpr Haugbolle& set_D( offset_t Ftype, int field_comp, int drv, Deriv_t<R,DGrid> f ) noexcept {
      _D[Ftype][field_comp][drv] = f;
      return *this;
    }
    constexpr Haugbolle& set_hh( offset_t Ftype, int k, HH_t<R> f ) noexcept {
      _hh[Ftype][k] = f;
      return *this;
    }
    constexpr Haugbolle& set_hh( offset_t Ftype, int field_comp, int drv, HH_t<R> f ) noexcept {
      return set_hh(!Ftype, 3 - field_comp - drv, f);
    }

    constexpr const auto& hh() const { return _hh; }
    constexpr const auto& D() const { return _D; }

    virtual Haugbolle* Clone() const { return new Haugbolle(*this); }

    virtual void operator() ( Field<R,3,DGrid>& E,
                              Field<R,3,DGrid>& B,
                              const Field<RJ,3,DGrid>& Jmesh,
                              const apt::Grid<R,DGrid>& grid,
                              const mpi::CartComm& cart,
                              int timestep,
                              R dt
                              ) const override;
  };

  template < typename R, int DGrid, typename RJ >
  struct HaugbolleBdry: public HaugbolleParams<HaugbolleBdry<R,DGrid,RJ>,R>,
                        public Action<R,DGrid,RJ> {
  private:
    R _preJ_factor {};
    int _num_iter = 4;

    Mesh<DGrid> _Dmesh;
    apt::array<apt::array< apt::array< std::vector<Deriv_t<R,DGrid>>,3 >,3>,2> _D{}; // derivatives D[Ftype][Fcomp][coordinate]
    apt::array< apt::array< std::vector< HH_t<R> >, 3 >, 2 > _hh {}; // hh[Ftype][Fcomp]

    apt::array<int,DGrid> _bdry{}; // -1 means lower, 0 means not a boundary, 1 means upper
    apt::array<bool,DGrid> _continuous_transverse_E{};

    int _idx_von_neumann ( apt::Index<DGrid> I ) const {
        for ( int i = 0; i < DGrid; ++i )
          I[i] = std::max( _Dmesh.range()[i].begin(), std::min( I[i], _Dmesh.range()[i].end() - 1 ) );
        return _Dmesh.linear_index(I);
      };

  public:
    HaugbolleBdry() = default;
    HaugbolleBdry(const HaugbolleBdry& ) = default;

    void init_from( const Haugbolle<R,DGrid,RJ>& haug ) {
      this->set_fourpi(haug.fourpi());
      this->set_number_iteration(haug.num_iter());
      this->set_implicit(haug.implicit());

      { // set up Dmesh
        apt::array<apt::Range,DGrid> r;
        for ( int i = 0; i < DGrid; ++i ) {
          r[i].begin() = apt::range::begin(*this,i);
          if ( _bdry[i] == 0 ) {
            r[i].end() = r[i].begin() + 1; // if not bdry, suppress the dimension
          } else {
            r[i].end() = apt::range::end(*this,i);
          }
        }
        _Dmesh = {r};
      }

      for ( int Ftype = 0; Ftype < 2; ++Ftype ) {
        for ( int c = 0; c < 3; ++c ) {
          _hh[Ftype][c].resize(_Dmesh.stride().back());
          for ( auto& x : _hh[Ftype][c] ) x = haug.hh()[Ftype][c];
        }

        for ( int c = 0; c < 3; ++c ) {
          for ( int i = 0; i < 3; ++i ) {
            if ( c == i ) continue; // only need transverse
            _D[Ftype][c][i].resize(_Dmesh.stride().back());
            for ( auto& x : _D[Ftype][c][i] ) x = haug.D()[Ftype][c][i];
          }
        }
      }
    }

    constexpr HaugbolleBdry& set_D( offset_t Ftype, int field_comp, int drv, Deriv_t<R,DGrid> f, const apt::Index<DGrid>& I ) noexcept {
      _D[Ftype][field_comp][drv][_idx_von_neumann(I)] = f;
      return *this;
    }

    constexpr HaugbolleBdry& set_hh( offset_t Ftype, int k, HH_t<R> f, const apt::Index<DGrid>& I ) noexcept {
      _hh[Ftype][k][_idx_von_neumann(I)] = f;
      return *this;
    }
    constexpr HaugbolleBdry& set_hh( offset_t Ftype, int field_comp, int drv, HH_t<R> f, const apt::Index<DGrid>& I ) noexcept {
      return set_hh(!Ftype, 3 - field_comp - drv, f, I);
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

    constexpr const auto& hh(offset_t Ftype, int k, const apt::Index<DGrid>& I) const
    { return _hh[Ftype][k][_idx_von_neumann(I)]; }
    constexpr const auto& D( offset_t Ftype, int field_comp, int drv, const apt::Index<DGrid>& I) const
    { return _D[Ftype][field_comp][drv][_idx_von_neumann(I)]; }
    constexpr const auto& boundary() const { return _bdry; }
    constexpr const auto& cont_trans_E() const { return _continuous_transverse_E; }

    virtual HaugbolleBdry* Clone() const { return new HaugbolleBdry(*this); }

    virtual void operator() ( Field<R,3,DGrid>& E,
                              Field<R,3,DGrid>& B,
                              const Field<RJ,3,DGrid>& Jmesh,
                              const apt::Grid<R,DGrid>& grid,
                              const mpi::CartComm& cart,
                              int timestep,
                              R dt
                              ) const override;
  };

}

#endif
