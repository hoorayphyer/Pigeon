#ifndef  _FIELD_CALCULUS_HPP_
#define  _FIELD_CALCULUS_HPP_

#include "field/field.hpp"

namespace field {
  template < knl::coordsys_t CS, typename T, int DGrid >
  struct calc {
  private:
    const knl::Grid<DGrid,T>& _grid;

  public:
    constexpr calc( const knl::Grid<DGrid,T>& grid ) noexcept : _grid(grid) {}

    void curl( Field<T,3,DGrid>& output, const Field<T,3,DGrid>& input );

    void div( Field<T,1,DGrid>& output, const Field<T,3,DGrid>& input );

    void lapl( Field<T,3,DGrid>& output, const Field<T,3,DGrid>& input );
  };

}

#endif
