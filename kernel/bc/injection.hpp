#ifndef _BC_INJECTION_HPP_
#define _BC_INJECTION_HPP_

namespace bc {
  template < int DGrid, typename Real >
  struct Injector {
    apt::Index<DGrid> Ib;
    apt::Index<DGrid> extent;

    Real v_th = 0.0;
    Real j_reg_x = 0.0;
    Real N_dot = 0.0;

    particle::species posion = particle::species::ion;
    particle::species negaon = particle::species::electron;

    Real (*omega_t) ( Real time ) = nullptr;

    template < template < typename > class Specs, typename RealJ >
    void operator() ( int timestep, Real dt, util::Rng<Real>& rng, particle::birthplace birth,
                      const mani::Grid<Real,DGrid>& grid, // local grid
                      const field::Field<Real, 3, DGrid>& E,
                      const field::Field<Real, 3, DGrid>& B,
                      const field::Field<RealJ, 3, DGrid>& J,
                      particle::map<particle::array<Real,Specs>>& particles
                      ) {
      using namespace particle;

      // the idea is that (timestep + 1) * N_inj * profile represents the accummulated number of injected pairs through the specified timestep. NOTE the timestep is shifted by one to reflect the actual times the inject is called, including this time.
      // NOTE _Jfield and Jmesh also differs in their underlying mesh guard cells
      auto profile_inj =
        [N_inj=N_dot*dt,timestep] ( Real theta ) noexcept {
          Real inj_num_base = N_inj * 0.5 * std::abs( std::sin( 2 * theta ) );
          return static_cast<int>( (timestep + 1) * inj_num_base ) - static_cast<int>( timestep * inj_num_base );
        };

      auto j_reg_inj =
        [factor = j_reg_x * mani::dV(grid) / (particle::properties.at(posion).charge_x) ]( Real J ) noexcept {
          return J * factor;
        };

      auto itr_po = std::back_inserter(particles[posion]);
      auto itr_ne = std::back_inserter(particles[negaon]);

      for ( auto I : apt::Block(extent) ) {
        I += Ib;
        apt::Vec<Real, Specs<Real>::Dim> q{};
        for ( int i = 0; i < DGrid; ++i )
          q[i] = grid[i].absc(I[i], 0.5);

        apt::Vec<Real,3> nB, J_v;
        // TODOL use interpolated value
        for ( int i = 0; i < 3; ++i ) {
          nB[i] = B[i](I);
          J_v[i] = J[i](I);
        }

        int num = std::max<Real>( profile_inj(q[1]), j_reg_inj(apt::abs(J_v)) );

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
          *(itr_ne++) = Particle<Real,Specs>( q_ptc, p_ptc, negaon, birth );
          *(itr_po++) = Particle<Real,Specs>( std::move(q_ptc), std::move(p_ptc), posion, birth );
        }
      }
    }
  };
}

#endif
