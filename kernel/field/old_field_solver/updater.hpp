#ifndef _OLDFIELDUPDATER_HPP_
#define _OLDFIELDUPDATER_HPP_

#include "field/field.hpp"
#include "field/action.hpp"
#include "apt/grid.hpp"

namespace mpi { struct CartComm; }

namespace field {
  template < typename R, int DGrid, typename RJ >
  struct OldSolve : public Action<R,DGrid,RJ> {
  private:
    R _fourpi {};
    const int _guard = 1; // just use this number
    int _magnetic_pole = 2; // 1 for mono-, 2 for di-
    int _surface_indent = 5;
    int _damp_indent = 43;
    R _damping_rate = 10.0;

    R (*_omega_t) (R);

    void Init( const mpi::CartComm& cart, const apt::Grid<R,DGrid>& local_grid) const;

  public:
    constexpr OldSolve& set_fourpi( R x ) noexcept {
      _fourpi = x;
      return *this;
    }

    constexpr OldSolve& set_magnetic_pole( int pole ) noexcept {
      _magnetic_pole = pole;
      return *this;
    }

    constexpr OldSolve& set_damping_rate( R rate ) noexcept {
      _damping_rate = rate;
      return *this;
    }

    constexpr OldSolve& set_surface_indent( int i ) noexcept {
      _surface_indent = i;
      return *this;
    }

    constexpr OldSolve& set_damp_indent( int i ) noexcept {
      _damp_indent = i;
      return *this;
    }

    constexpr OldSolve& set_omega_t( R(*f)(R) ) noexcept {
      _omega_t = f;
      return *this;
    }

    constexpr int guard() const noexcept { return _guard; }

    virtual OldSolve* Clone() const { return new OldSolve(*this); }

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
