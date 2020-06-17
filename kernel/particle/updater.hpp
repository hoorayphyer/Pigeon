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

#include "particle/migration.hpp"

namespace particle {

  template < int DGrid,
             typename R,
             template < typename > class S,
             typename ShapeF,
             typename RJ >
  class Migrator : public Action<DGrid,R,S,RJ> {
  private:
    apt::Grid< R, DGrid > _supergrid;

  public:
    Migrator& set_supergrid( const apt::Grid< R, DGrid >& supergrid ) noexcept { _supergrid = supergrid; return *this; }
    Migrator* Clone() const override { return new Migrator(*this); }

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
      // bulk range = [lb, ub)
      constexpr auto migrtrit =
        []( auto q, auto lb, auto ub ) noexcept {
          return ( q >= lb ) + ( q >= ub );
        };

      for ( auto sp : particles ) {
        for ( auto ptc : particles[sp] ) { // TODOL semantics
          if ( !ptc.is(flag::exist) ) continue;
          int mig_dir{};
          for ( int i = 0; i < DGrid; ++i ) {
            mig_dir += migrtrit( ptc.q(i), grid[i].lower(), grid[i].upper() ) * apt::pow3(i);
          }

          if ( mig_dir != ( apt::pow3(DGrid) - 1 ) / 2 ) {
            ptc.template set<migrcode,DGrid>(mig_dir);
            new_ptc_buf->emplace_back(std::move(ptc));
          }
        }
      }

      migrate( *new_ptc_buf, ens->cart_topos, ens->inter, timestep );

      // NOTE adjust particle positions in the ring topology, regardless how many cpus there are on that ring.
      {
        bool has_periodic = false;
        for ( int i = 0; i < DGrid; ++i ) {
          has_periodic = has_periodic || ens->cart_topos[i].periodic();
        }
        if ( has_periodic ) {
          for ( auto& ptc : *new_ptc_buf ) {
            if ( !ptc.is(flag::exist) ) continue;

            for ( int i = 0; i < DGrid; ++i ) {
              if ( !ens->cart_topos[i].periodic() ) continue;
              int idx = static_cast<int>( ( ptc.q(i) - _supergrid[i].lower() ) / _supergrid[i].delta() + 0.5 );
              if ( idx >= 0 ) idx /= _supergrid[i].dim();
              else idx = - ( (-idx) / _supergrid[i].dim() + 1 );

              ptc.q(i) -= idx * _supergrid[i].dim() * _supergrid[i].delta();
            }
          }
        }
      }

      for ( auto&& ptc : *new_ptc_buf ) {
        if ( !ptc.is(flag::exist) ) continue;
        auto sp = ptc.template get<species>();
#ifdef PIC_DEBUG
        // // check if the received ptc trully resides in this ensemble.
        // apt::array<int,DGrid> mig_co;
        // bool is_OK = true;
        // for ( int i = 0; i < DGrid; ++i ) {
        //   mig_co[i] = migrate_code( ptc.q(i), _borders[i][LFT], _borders[i][RGT] );
        //   if ( mig_co[i] != 1 && !_ens_opt->is_at_boundary(i)[(mig_co[i] != 0)] ) //NOTE need to consider boundaries
        //     is_OK = false;
        // }
        // if ( !is_OK ) {
        //   lgr::file << "ts=" << debug::timestep << ", wr=" << debug::world_rank << ", el=" << debug::ens_label << std::endl;
        //   lgr::file << "Received across-ensemble particles! q = " << ptc.q() << ", p = " << ptc.p() << std::endl;
        //   lgr::file << "  mig_dir on new ensemble  = " << mig_co;
        //   // get old mig_co
        //   for ( int i = 0; i < DGrid; ++i ) {
        //     mig_co[i] = ( migrInt<DGrid>(ptc) % apt::pow3(i+1) ) / apt::pow3(i);
        //   }
        //   lgr::file << ", mig_dir on old ensemble = " << mig_co << std::endl;
        //   debug::throw_error("Received across-ensemble particles!");
        // }
#endif
        ptc.template reset<migrcode>();
        particles[sp].push_back( std::move(ptc) );
      }
      new_ptc_buf->resize(0);
    }
  };

}

#endif
