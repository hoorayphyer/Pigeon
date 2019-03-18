#ifndef _PARTICLE_GRID_BASED_INJECTION_HPP_
#define _PARTICLE_GRID_BASED_INJECTION_HPP_

#include "apt/block.hpp"
#include "apt/apply.hpp"

namespace particle {
  template < typename InjectionNumberPolicy,
             typename NegatonInjectPolicy,
             typename PositonInjectPolicy >
  struct GridBasedInjector {
  private:
    const InjectionNumberPolicy& _num_inj;
    const NegatonInjectPolicy& _inj_negaton;
    const NegatonInjectPolicy& _inj_positon;

  public:
    constexpr GridBasedInjector( const InjectionNumberPolicy& num_inj,
                                 const NegatonInjectPolicy& inj_negaton,
                                 const NegatonInjectPolicy& inj_positon,
                                 ) noexcept
      : _num_inj(num_inj), _inj_negaton(inj_negaton), _inj_positon(inj_positon) {}

    template < typename BackInsertIter, typename Grid >
    void inject( BackInsertIter itr_negaton,
                 BackInsertIter itr_positon,
                 const Grid& grid ) {
      // TODO over-injection protection? Do it completely in the coordinate space
      // TODO actual injection is reduced by a factor of ensemble size, so is coordinate current

      typename T = Grid::element_t;
      constexpr int D = Grid::NDim;
      apt::Index<D> extent;

      apt::foreach<0, D>
        ( [](auto& e, const auto& g) noexcept
          { e = g.dim(); }, extent, grid );

      apt::array<T, D> loc{};

      for ( const auto& I : apt::block(extent) ) {
        apt::foreach<0,D>
          ( [](auto& l, const auto& g, auto i ) noexcept
            {l = g.absc(i, 0.5);}, loc, grid, I );

        unsigned int num = apt::apply( _num_inj, I );

        while ( num-- ) {
          *(itr_negaton++) = _inj_negaton(loc);
          *(itr_positon++) = _inj_positon(loc);
        }
      }
    }

  };

}


#endif
