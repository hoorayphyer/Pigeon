#ifndef _FIELD_SHAPE_INTERPLAY_HPP_
#define _FIELD_SHAPE_INTERPLAY_HPP_

#include "field/field.hpp"
#include "apt/vec.hpp"

namespace mpi { struct Comm; struct CartComm; }

namespace field {
  template < typename T, int DField, int DGrid >
  class dJ_Field{
  private:
    Field<T,DField,DGrid> _data;

  public:
    dJ_Field( const Mesh<T,DGrid>& mesh );

    template < typename Vec_q0, typename Vec_q1, typename ShapeF >
    void deposit ( T charge_over_dt,
                   const apt::VecExpression<Vec_q0>& q0_abs,
                   const apt::VecExpression<Vec_q1>& q1_abs,
                   const ShapeF& shapef );

    void reduce( int chief, const mpi::Comm& intra );

    const Field<T,DField,DGrid>& integrate( const mpi::CartComm& cart );

    inline void reset() {
      apt::foreach<0, DField>
        ( []( auto comp ) { // returns a proxy
            for ( auto& elm : comp.data() ) elm = 0.0;
          }, _data );
    }

  };

}

namespace field {

  template < typename T, int DField, int DGrid, typename LocType, typename ShapeF >
  apt::Vec<T, DField> interpolate ( const Field<T,DField,DGrid>& field,
                                    const LocType& q_abs,
                                    const ShapeF& shapef );

  // the opposite of interpolate. `var` is +=ed to field
  template < typename T, int DField, int DGrid, typename LocType, typename ShapeF >
  void deposit ( Field<T,DField,DGrid>& field,
                 apt::Vec<T, DField> var,
                 const LocType& q_abs,
                 const ShapeF& shapef );
}
#endif
