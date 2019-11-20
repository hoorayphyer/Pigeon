#ifndef  _PARTICLE_UPDATER_HPP_
#define  _PARTICLE_UPDATER_HPP_

#include "particle/action.hpp"
#include "particle/species_predef.hpp"
#include "particle/forces.hpp"
#include "particle/scattering.hpp"
#include "particle/load_type.hpp"

namespace dye {
  template < int > struct Ensemble;
}

namespace particle {
  template < int DGrid,
             typename R,
             template < typename > class S,
             typename ShapeF,
             typename RJ >
  class Updater : public Action<DGrid,R,S,RJ> {
  private:
    apt::array<R, S<R>::Dim> (*_update_q)( typename array<R,S>::particle_type::vec_type& x, typename array<R,S>::particle_type::vec_type& p, R dt, bool is_massive );

  public:
    Updater* Clone() const override { return new Updater(*this); }
    Updater& set_update_q( apt::array<R, S<R>::Dim> (*update_q)( typename array<R,S>::particle_type::vec_type&, typename array<R,S>::particle_type::vec_type&, R, bool ) ) { _update_q = update_q; return *this; }

    void operator() ( map<array<R,S>>& particles,
                      field::Field<RJ,3,DGrid>& J,
                      std::vector<Particle<R,S>>* new_ptc_buf,
                      const map<Properties>& properties,
                      const field::Field<R,3,DGrid>& E,
                      const field::Field<R,3,DGrid>& B,
                      const apt::Grid< R, DGrid >& grid,
                      const dye::Ensemble<DGrid>* ,
                      R dt, int timestep, util::Rng<R>& rng
                      ) override;
  };

}

#endif
