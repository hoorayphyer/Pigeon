#ifndef  _PARTICLE_ANNIHILATION_HPP_
#define  _PARTICLE_ANNIHILATION_HPP_

#include "particle/action.hpp"
#include "dye/ensemble.hpp"
#include <cassert>

namespace particle {
  template < int DGrid, typename R, template < typename > class S, typename ShapeF, typename RJ >
  void annihilate( array<R,S>& el, array<R,S>& po,
                   field::Field<RJ,3,DGrid>& J,
                   R charge_el, R charge_po,
                   const apt::Grid< R, DGrid >& grid,
                   const mpi::Comm& intra,
                   R dt, const ShapeF&,
                   R(*policy)(R num_electron_in_a_cell, R num_positron_in_a_cell) );

  template < int DGrid,
             typename R,
             template < typename > class S,
             typename ShapeF,
             typename RJ >
  class Annihilator : public Action<DGrid,R,S,RJ> {
  private:
    R(*_policy)(R num_electron_in_a_cell, R num_positron_in_a_cell) = nullptr;

  public:
    Annihilator* Clone() const override { return new Annihilator(*this); }

    Annihilator& set_policy( R(*policy)(R,R) ) noexcept { _policy = policy; return *this; }

    void operator() ( map<array<R,S>>& particles,
                      field::Field<RJ,3,DGrid>& J,
                      std::vector<Particle<R,S>>* new_ptc_buf,
                      const map<Properties>& properties,
                      const field::Field<R,3,DGrid>& E,
                      const field::Field<R,3,DGrid>& B,
                      const apt::Grid< R, DGrid >& grid,
                      const dye::Ensemble<DGrid>* ens,
                      R dt, int timestep, util::Rng<R>& rng
                      ) override {
      if ( apt::range::is_empty(*this) ) return;
      assert( ens != nullptr );
      assert( _policy != nullptr );

      annihilate(particles[species::electron],particles[species::positron],J,
                 properties[species::electron].charge_x,properties[species::positron].charge_x,
                 grid, ens->intra, dt, ShapeF(), _policy);
    }

  };
}

#endif
