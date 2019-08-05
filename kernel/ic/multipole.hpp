#ifndef _IC_MULTIPOLE_HPP_
#define _IC_MULTIPOLE_HPP_

namespace ic {
  template < int DGrid >
  struct MagneticDipole {
    apt::Index<DGrid> Ib;
    apt::Index<DGrid> extent;

    template < typename T >
    T B_r ( T logr, T theta ) noexcept {
      return 2.0 * std::cos(theta) * std::exp(-3.0 * logr);
    };

    template < typename T >
    T B_th ( T logr, T theta ) noexcept {
      return std::sin(theta) * std::exp(-3.0 * logr);
    };

    template < typename Real, typename RealJ, template < typename > class Specs >
    void operator() ( const mani::Grid<Real,DGrid>& grid,
                      field::Field<Real, 3, DGrid>& E,
                      field::Field<Real, 3, DGrid>& B,
                      field::Field<RealJ, 3, DGrid>& J,
                      particle::map<particle::array<Real,Specs>>&
                      ) {
      for ( auto I : apt::Block(extent) ) {
        I += Ib;
        B[0](I) = B_r( grid[0].absc(I[0], B[0].offset()[0]), grid[1].absc(I[1], B[0].offset()[1]) );
        B[1](I) = B_th( grid[0].absc(I[0], B[1].offset()[0]), grid[1].absc(I[1], B[1].offset()[1]) );
      }
    }

    constexpr int initial_timestep() noexcept { return 0; }
  };
}

#endif
