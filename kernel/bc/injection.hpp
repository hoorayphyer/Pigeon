#ifndef _BC_INJECTION_HPP_
#define _BC_INJECTION_HPP_

namespace bc {
  template < int DGrid, typename Real >
  struct Injector {
  private:
    field::Field<Real,1,DGrid> _count_n;
    field::Field<Real,1,DGrid> _count_p;

  public:
    apt::Index<DGrid> Ib;
    apt::Index<DGrid> extent;

    Real v_th = 0.0;
    Real N_atm = 0.0;

    particle::species posion = particle::species::ion;
    particle::species negaon = particle::species::electron;

    Real (*omega_t) ( Real time ) = nullptr;

    template < template < typename > class Specs, typename RealJ >
    void operator() ( int timestep, Real dt, util::Rng<Real>& rng,
                      Real wdt_pic,
                      const dye::Ensemble<DGrid>& ens,
                      const mani::Grid<Real,DGrid>& grid, // local grid
                      const field::Field<Real, 3, DGrid>& E,
                      const field::Field<Real, 3, DGrid>& B,
                      const field::Field<RealJ, 3, DGrid>& J,
                      particle::map<particle::array<Real,Specs>>& particles
                      ) {
      using namespace particle;

      _count_n.resize( {extent,0} );
      _count_p.resize( {extent,0} );

      apt::array<Real,DGrid> lb;
      apt::array<Real,DGrid> ub;
      for ( int i = 0; i < DGrid; ++i ) {
        lb[i] = grid[i].absc( Ib[i], 0.0 );
        ub[i] = grid[i].absc( Ib[i] + extent[i], 0.0 );
      }
      auto is_in = [&lb,&ub]( const auto& q ) {
                     for ( int i = 0; i < DGrid; ++i ) {
                       if ( q[i] < lb[i] || q[i] >= ub[i] ) return false;
                     }
                     return true;
                   };

      auto f_count
        = [&lb,&ub,&grid,is_in]( auto& count, const auto& ptcs) {
            count.reset();
            for ( const auto& x : ptcs ) {
              if ( !x.is(flag::exist) || !is_in(x.q()) ) continue;
              apt::Index<DGrid> idx;
              for ( int i = 0; i < DGrid; ++i )
                idx[i] = ( x.q()[i] - lb[i] ) / grid[i].delta();
              count[0](idx) += x.frac(); // add by fraction
            }
          };

      f_count( _count_n, particles[negaon] );
      f_count( _count_p, particles[posion] );

      { // parallelizae TODO optimize
        int rank_inj = timestep % ens.size();
        ens.intra.template reduce<true>(mpi::by::SUM, rank_inj, _count_n[0].data().data(), _count_n[0].data().size() );
        ens.intra.template reduce<true>(mpi::by::SUM, rank_inj, _count_p[0].data().data(), _count_p[0].data().size() );
        if ( ens.intra.rank() != rank_inj ) return;
      }

      auto itr_po = std::back_inserter(particles[posion]);
      auto itr_ne = std::back_inserter(particles[negaon]);

      for ( auto I : apt::Block(extent) ) {
        auto N_pairs = std::min( _count_n[0](I), _count_p[0](I) );
        I += Ib;
        apt::Vec<Real, Specs<Real>::Dim> q{};
        for ( int i = 0; i < DGrid; ++i )
          q[i] = grid[i].absc(I[i], 0.5);

        apt::Vec<Real,3> nB;
        { // make nB centered in the cell
          const auto& m = B.mesh();
          auto li = m.linearized_index_of_whole_mesh(I);
          if constexpr (DGrid == 2) {
              nB[0] = 0.5 * ( B[0][li] + B[0][li + m.stride(1)] );
              nB[1] = 0.5 * ( B[1][li] + B[1][li + m.stride(0)] );
              nB[2] = 0.25 * ( B[2][li] + B[2][li + m.stride(0)] + B[2][li + m.stride(1)] + B[2][li + m.stride(0) + m.stride(1)] );
                } else if (DGrid == 3){
            nB[0] = 0.25 * ( B[0][li] + B[0][li + m.stride(1)] + B[0][li + m.stride(2)] + B[0][li + m.stride(1) + m.stride(2)] );
            nB[1] = 0.25 * ( B[1][li] + B[1][li + m.stride(2)] + B[1][li + m.stride(0)] + B[1][li + m.stride(2) + m.stride(0)] );
            nB[2] = 0.25 * ( B[2][li] + B[2][li + m.stride(0)] + B[2][li + m.stride(1)] + B[2][li + m.stride(0) + m.stride(1)] );
          }
          if ( apt::abs(nB) == 0.0 ) nB = {1.0, 0.0, 0.0}; // use radial direction as default
          else nB /= apt::abs(nB);
        }

        apt::Vec<Real, Specs<Real>::Dim> p{};
        p[2] = omega_t( timestep * dt ) * std::exp(q[0]) * std::sin(q[1]); // corotating

        // replenish
        Real quota = N_atm * std::sin(q[1]) - N_pairs;
        while ( quota > 0 ) {
          auto q_ptc = q;
          Real frac = std::min( (Real)1.0, quota );
          quota -= (Real)1.0;

          for ( int i = 0; i < DGrid; ++i )
            q_ptc[i] += grid[i].delta() * rng.uniform(-0.5, 0.5);
          auto p_ptc = p;
          p_ptc += nB * rng.gaussian( 0.0, v_th );
          *(itr_ne++) = Particle<Real,Specs>( q_ptc, p_ptc, frac, negaon, birthplace(ens.label()) );
          *(itr_po++) = Particle<Real,Specs>( std::move(q_ptc), std::move(p_ptc), frac, posion, birthplace(ens.label()) );
        }
      }
    }
  };
}

#endif
