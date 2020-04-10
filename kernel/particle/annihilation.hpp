#ifndef  _PARTICLE_ANNIHILATION_HPP_
#define  _PARTICLE_ANNIHILATION_HPP_
#include "particle/action.hpp"

namespace particle {
  template < int DGrid,
             typename R,
             template < typename > class S,
             typename ShapeF,
             typename RJ >
  class Annihilator : public Action<DGrid,R,S,RJ> {
  private:

  public:
    Annihilator* Clone() const override { return new Annihilator(*this); }

    void operator() ( map<array<R,S>>& particles,
                      field::Field<RJ,3,DGrid>& J,
                      std::vector<Particle<R,S>>* new_ptc_buf,
                      const map<Properties>& properties,
                      const field::Field<R,3,DGrid>& E,
                      const field::Field<R,3,DGrid>& B,
                      const apt::Grid< R, DGrid >& grid,
                      const dye::Ensemble<DGrid>* ens,
                      R dt, int timestep, util::Rng<R>& rng
                      ) override;
  };
}

#endif
