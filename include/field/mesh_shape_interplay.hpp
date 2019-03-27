#ifndef _FIELD_MESH_SHAPE_INTERPLAY_HPP_
#define _FIELD_MESH_SHAPE_INTERPLAY_HPP_

#include "field/field.hpp"
#include "apt/vec.hpp"

namespace mpi { struct Comm; struct CartComm; }

namespace field {
  // NOTE in Standard_dJ_Field, the standard grid( or bulk of the mesh) has lower bound index = 0 and spacing = 1. "Standard" is borrowed from "standard normal distribution". This way, there is no need of grid information. NOTE though the dJ field is still defined in the MIDWAY.
  template < typename T, int DField, int DGrid >
  class Standard_dJ_Field{
  private:
    Field<T,DField,DGrid> _data;

  public:
    Standard_dJ_Field( const Mesh<DGrid>& mesh );

    // TODO double check if calculation precision is promoted besides assignment
    template < typename Vec_q0, typename Vec_q1, typename ShapeF >
    void deposit ( T charge_over_dt,
                   const apt::VecExpression<Vec_q0>& q0_std,
                   const apt::VecExpression<Vec_q1>& q1_std,
                   const ShapeF& shapef );

    void reduce( int chief, const mpi::Comm& intra );

    Field<T,DField,DGrid>& integrate( const mpi::CartComm& cart );

    inline void reset() {
      apt::foreach<0, DField>
        ( []( auto comp ) { // returns a proxy
            for ( auto& elm : comp.data() ) elm = 0.0;
          }, _data );
    }

  };

}

namespace field {

  // NOTE q_std refers to the same "standard" as above
  template < typename T, int DField, int DGrid, typename ShapeF >
  apt::Vec<T, DField> interpolate ( const Field<T,DField,DGrid>& field,
                                    const apt::array<T,DGrid>& q_std,
                                    const ShapeF& shapef );

  // the opposite of interpolate. `var` is +=ed to field
  template < typename T, int DField, int DGrid, typename ShapeF >
  void deposit ( Field<T,DField,DGrid>& field,
                 apt::Vec<T, DField> var,
                 const apt::array<T,DGrid>& q_std,
                 const ShapeF& shapef );
}
#endif
