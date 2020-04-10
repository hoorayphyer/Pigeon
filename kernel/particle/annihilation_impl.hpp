#include "particle/annihilation.hpp"
#include <cassert>
#include "msh/mesh_shape_interplay.hpp"
#include "mpipp/mpi++.hpp"
#include "dye/ensemble.hpp"

namespace particle {
  template < int DGrid, typename R, template < typename > class S, typename ShapeF, typename RJ >
  void annihilate( array<R,S>& el, array<R,S>& po,
                   field::Field<RJ,3,DGrid>& J,
                   R charge_el, R charge_po,
                   const apt::Grid< R, DGrid >& grid,
                   const mpi::Comm& intra,
                   R dt, const ShapeF& ) {
    const auto& mesh = J.mesh();
    const auto size = mesh.linear_size();
    std::vector<R> buf (2 * size, 0);

    auto get_idx =
      [&] ( const auto& ptc ) {
        apt::Index<DGrid> I;
        for ( int i = 0; i < DGrid; ++i )
          I[i] = ( ptc.q(i) - grid[i].lower() ) / grid[i].delta();
        return mesh.linear_index(I);
      };

    {
      std::vector<R> annih (size, 0);

      auto f =
        [&] (const auto& ptcs, auto* ptr) {
          for ( const auto& x : ptcs) {
            if ( not x.is(flag::exist) ) continue;
            ptr[get_idx(x)] += x.frac();
          }
        };
      f(el, buf.data());
      f(po, buf.data() + size);

      if ( intra.size() > 1 )
        intra.inscan_inplace(mpi::by::SUM, buf.data(), buf.size());

      if ( intra.rank() == intra.size() - 1 ) {
        for ( auto i = decltype(size){0}; i < size; ++i ) {
          annih[i] = std::min(buf[i], buf[size+i]);
        }
      }

      if (intra.size() > 1)
        intra.broadcast( intra.size()-1, annih.data(), annih.size() );

      for ( auto i = decltype(size){0}; i < size; ++i )
        buf[i] -= annih[i];

      for ( auto i = decltype(size){0}; i < size; ++i )
        buf[size+i] -= annih[i];
    }
    {
      constexpr auto shapef = ShapeF();

      auto f =
        [&] (auto& ptcs, auto* ptr, auto charge_over_dt) {
          for ( auto x : ptcs) { // TODOL semantics
            if ( not x.is(flag::exist) ) continue;
            auto& exempt = ptr[get_idx(x)];
            if ( exempt > x.frac() ) {
              exempt -= x.frac();
            } else {
              bool is_partial = (exempt > static_cast<R>(0));
              if ( is_partial ) x.frac() -= exempt;

              auto q0_std = msh::to_standard( grid, x.q() );
              decltype(q0_std) q1_std;
              for ( int i = 0; i < DGrid; ++i ) {// NOTE DGrid is used here, not DPtc
                q1_std[i] = std::floor(q0_std[i]) + static_cast<R>(0.5);
              }
              for ( int i = DGrid; i < S<R>::Dim; ++i ) {// NOTE extra dimensions, no current is generated
                q1_std[i] = q0_std[i];
              }
              msh::deposit( J, x.frac() * charge_over_dt, shapef, q0_std, q1_std );

              if ( is_partial ) {
                x.frac() = exempt;
                exempt = static_cast<R>(0);
              } else {
                x.reset(flag::exist);
              }
            }

          }
        };

      f(el,buf.data(),charge_el / dt);
      f(po,buf.data()+size,charge_po / dt);
    }
  }


  template < int DGrid, typename R, template < typename > class S, typename ShapeF, typename RJ >
  void Annihilator<DGrid,R,S,ShapeF,RJ>::operator() ( map<array<R,S>>& particles,
                                                      field::Field<RJ,3,DGrid>& J,
                                                      std::vector<Particle<R,S>>* new_ptc_buf,
                                                      const map<Properties>& properties,
                                                      const field::Field<R,3,DGrid>& E,
                                                      const field::Field<R,3,DGrid>& B,
                                                      const apt::Grid< R, DGrid >& grid,
                                                      const dye::Ensemble<DGrid>* ens,
                                                      R dt, int timestep, util::Rng<R>& rng
                                                      ) {
    assert( ens != nullptr );

    annihilate(particles[species::electron],particles[species::positron],J,
               properties[species::electron].charge_x,properties[species::positron].charge_x,
               grid, ens->intra, dt, ShapeF());
  }
}
