#ifndef _BC_AXISSYMMETRIC_HPP_
#define _BC_AXISSYMMETRIC_HPP_

namespace bc {
  // template < int DGrid,
  //            typename Real,
  //            template < typename > class Specs,
  //            typename RealJ >
  // struct FieldBC_Rotating_Conductor {
  // private:
  //   const mani::Grid<Real,DGrid>& _grid;
  //   field::Field<Real, 3, DGrid>& _Efield;
  //   field::Field<Real, 3, DGrid>& _Bfield;

  //   apt::Index<DGrid> _Ib;
  //   apt::Index<DGrid> _extent;

  //   const Real _mu0 = pic::mu0;

  //   Real B_r_over_mu0 ( Real logr, Real theta ) noexcept {
  //     return 2.0 * std::cos(theta) * std::exp(-3.0 * logr);
  //   };

  //   Real B_th_over_mu0 ( Real logr, Real theta ) noexcept {
  //     return std::sin(theta) * std::exp(-3.0 * logr);
  //   };

  //   // find out E by E = -( Omega x r ) x B
  //   Real E_r_over_mu0omega ( Real logr, Real theta ) noexcept {
  //     auto sin_t = std::sin(theta);
  //     return std::exp( -2.0 * logr ) * sin_t * sin_t;
  //   };

  //   Real E_th_over_mu0omega ( Real logr, Real theta ) noexcept {
  //     return - std::exp( -2.0 * logr ) * std::sin( 2.0 * theta );
  //   };

  // public:
  //   FieldBC_Rotating_Conductor ( const mani::Grid<Real,DGrid>& localgrid,
  //                                field::Field<Real, 3, DGrid>& Efield,
  //                                field::Field<Real, 3, DGrid>& Bfield,
  //                                const field::Field<RealJ, 3, DGrid>& Jfield, // J is Jmesh on a replica
  //                                const particle::map<particle::array<Real,Specs>>& particles )
  //     : _grid(localgrid), _Efield(Efield), _Bfield(Bfield) {
  //     _Ib[0] = 0;
  //     _extent[0] = pic::ofs::indent[0];
  //     apt::tie(_Ib[1], _extent[1]) = gtl( {0.0, PI}, localgrid[1] );
  //   }

  //   void operator() ( int timestep, Real dt ) {
  //     const auto omega = pic::omega_spinup( timestep * dt );
  //     for ( auto I : apt::Block(_extent) ) {
  //       I += _Ib;
  //       // TODO deal with discontinuous and continuous variables
  //       // TODO piecing together different patches
  //       _Bfield[0](I) = _mu0 * B_r_over_mu0( _grid[0].absc(I[0], _Bfield[0].offset()[0]), _grid[1].absc(I[1], _Bfield[0].offset()[1]) );
  //       _Bfield[1](I) = _mu0 * B_th_over_mu0( _grid[0].absc(I[0], _Bfield[1].offset()[0]), _grid[1].absc(I[1], _Bfield[1].offset()[1]) );
  //       _Bfield[2](I) = 0.0;

  //       _Efield[0](I) = _mu0 * omega * E_r_over_mu0omega( _grid[0].absc(I[0], _Efield[0].offset()[0]), _grid[1].absc(I[1], _Efield[0].offset()[1]) );
  //       _Efield[1](I) = _mu0 * omega * E_th_over_mu0omega( _grid[0].absc(I[0], _Efield[1].offset()[0]), _grid[1].absc(I[1], _Efield[1].offset()[1]) );
  //       _Efield[2](I) = 0.0;
  //     }
  //   }
  // };

}

#endif
