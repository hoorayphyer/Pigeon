#include "particle/annihilation.hpp"
#include "msh/mesh_shape_interplay.hpp"
#include "mpipp/mpi++.hpp"
#include <cmath>

namespace particle {
  template < int DGrid, typename R, template < typename > class S, typename ShapeF, typename RJ >
  void annihilate( array<R,S>& el, array<R,S>& po,
                   field::Field<RJ,3,DGrid>& J,
                   R charge_el, R charge_po,
                   const apt::Grid< R, DGrid >& grid,
                   const mpi::Comm& intra,
                   R dt, const ShapeF&,
                   R(*policy)(R num_electron_in_a_cell, R num_positron_in_a_cell),
                   const apt::array<apt::array<R,2>,DGrid>& bounds ) {
    const auto& mesh = J.mesh();
    { // check if bounds is applicable
      apt::array<apt::Range,DGrid> range;
      for ( int i = 0; i < DGrid; ++i ) {
        range[i].begin() = ( bounds[i][0] - grid[i].lower() ) / grid[i].delta();
        range[i].end() = ( bounds[i][1] - grid[i].lower() ) / grid[i].delta();

        // NOTE far_begin/end are used
        range[i].begin() = std::max( range[i].begin(), mesh.range(i).far_begin() );
        range[i].end() = std::min( range[i].end(), mesh.range(i).far_end() );
      }

      if ( apt::range::is_empty(range) ) return;
    }

    const auto size = mesh.linear_size();
    std::vector<R> buf (2 * size, 0);

    auto get_idx =
      [&] ( const auto& ptc ) {
        apt::Index<DGrid> I;
        for ( int i = 0; i < DGrid; ++i )
          I[i] = ( ptc.q(i) - grid[i].lower() ) / grid[i].delta();
        return mesh.linear_index(I);
      };

    auto eligible =
      [&bounds] ( const auto& q ) noexcept {
        for ( int i = 0; i < DGrid; ++i ) {
          if ( q[i] < bounds[i][0] or q[i] >= bounds[i][1] )
            return false;
        }
        return true;
      };

    {
      std::vector<R> annih (size, 0);

      auto f =
        [&] (const auto& ptcs, auto* ptr) {
          for ( const auto& x : ptcs) {
            if ( not x.is(flag::exist) or not eligible(x.q()) ) continue;
            ptr[get_idx(x)] += x.frac();
          }
        };
      f(el, buf.data());
      f(po, buf.data() + size);

      if ( intra.size() > 1 )
        intra.inscan_inplace(mpi::by::SUM, buf.data(), buf.size());

      if ( intra.rank() == intra.size() - 1 ) {
        for ( auto i = decltype(size){0}; i < size; ++i ) {
          annih[i] = policy( buf[i], buf[size+i]  );
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
            if ( not x.is(flag::exist) or not eligible(x.q()) ) continue;
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

}
