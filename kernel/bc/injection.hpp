#ifndef _BC_INJECTION_HPP_
#define _BC_INJECTION_HPP_

namespace bc {
  template < int DGrid, typename Real >
  struct Injector {
  private:
    field::Field<int,1,DGrid> _count_n;
    field::Field<int,1,DGrid> _count_p;

  public:
    apt::Index<DGrid> Ib;
    apt::Index<DGrid> extent;

    Real v_th = 0.0;
    Real j_reg_x = 0.0;
    Real N_dot = 0.0;

    particle::species posion = particle::species::ion;
    particle::species negaon = particle::species::electron;

    Real (*omega_t) ( Real time ) = nullptr;

    template < template < typename > class Specs, typename RealJ >
    void operator() ( int timestep, Real dt, util::Rng<Real>& rng,
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
              ++count[0](idx);
            }
          };

      f_count( _count_n, particles[negaon] );
      f_count( _count_p, particles[posion] );

      // the idea is that (timestep + 1) * N_inj * profile represents the accummulated number of injected pairs through the specified timestep. NOTE the timestep is shifted by one to reflect the actual times the inject is called, including this time.
      // NOTE _Jfield and Jmesh also differs in their underlying mesh guard cells
      auto profile_inj =
        [N_inj=N_dot*dt,timestep, &ens] ( Real theta ) noexcept {
          Real inj_num_base = N_inj * 0.5 * std::abs( std::sin( 2 * theta ) );
          auto N_tot = static_cast<int>( (timestep + 1) * inj_num_base ) - static_cast<int>( timestep * inj_num_base );
          int r = ens.intra.rank();
          int size = ens.size();
          return N_tot / size + ( ((r + timestep) % size) < (N_tot % size) );
        };

      auto j_reg_inj =
        [factor = j_reg_x * mani::dV(grid) / (particle::properties.at(posion).charge_x) ]( Real J ) noexcept {
          return J * factor;
        };

      auto itr_po = std::back_inserter(particles[posion]);
      auto itr_ne = std::back_inserter(particles[negaon]);

      for ( auto I : apt::Block(extent) ) {
        int N_pairs = std::min( _count_n[0](I), _count_p[0](I) );
        I += Ib;
        apt::Vec<Real, Specs<Real>::Dim> q{};
        for ( int i = 0; i < DGrid; ++i )
          q[i] = grid[i].absc(I[i], 0.5);

        apt::Vec<Real,3> nB, J_v;
        for ( int i = 0; i < 3; ++i ) {
          nB[i] = B[i](I);
          J_v[i] = J[i](I);
        }

        int num = std::max<Real>( profile_inj(q[1]), j_reg_inj(apt::abs(J_v)) - N_pairs );

        // find n_B
        nB /= apt::abs(nB);

        apt::Vec<Real, Specs<Real>::Dim> p{};
        p[2] = omega_t( timestep * dt ) * std::exp(q[0]) * std::sin(q[1]); // corotating

        for ( int n = 0; n < num; ++n ) {
          auto q_ptc = q;
          for ( int i = 0; i < DGrid; ++i )
            q_ptc[i] += grid[i].delta() * rng.uniform(-0.5, 0.5);
          auto p_ptc = p;
          p_ptc += nB * rng.gaussian( 0.0, v_th );
          *(itr_ne++) = Particle<Real,Specs>( q_ptc, p_ptc, negaon, birthplace(ens.label()) );
          *(itr_po++) = Particle<Real,Specs>( std::move(q_ptc), std::move(p_ptc), posion, birthplace(ens.label()) );
        }
      }
    }
  };
}

#endif
