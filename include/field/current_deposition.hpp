#ifndef _FIELD_CURRENT_DEPOSITION_HPP_
#define _FIELD_CURRENT_DEPOSITION_HPP_

#include "field/field.hpp"
#include "apt/array.hpp"

namespace mpi { struct Comm; struct CartComm; }

namespace field {
  // NOTE The standard grid of the bulk of the mesh is a rescaled and shifted version of the actual grid such that the spacing = 1 and the first cell has index = 0. This way, there is no need of grid information. "Standard" is borrowed from "standard normal distribution".
  // NOTE the dJ field is still defined in the MIDWAY.
  template < typename T, int DField, int DGrid, typename ShapeF >
  class Standard_dJ_Field{
  private:
    Field<T,DField,DGrid> _data;

  public:
    Standard_dJ_Field( apt::Index<DGrid> bulk_extent, const ShapeF& );

    template < typename U >
    void deposit ( U charge_over_dt,
                   const apt::array<U,DField>& q0_std, // NOTE its DField, not DGrid
                   const apt::array<U,DField>& q1_std );

    void reduce( int chief, const mpi::Comm& intra );

    Field<T,DField,DGrid>& integrate( const mpi::CartComm& cart );

    inline void reset() {
      apt::foreach<0, DField>
        ( []( auto comp ) { // TODOL semantics
            for ( auto& elm : comp.data() ) elm = 0.0;
          }, _data );
    }

    inline const auto& mesh() const noexcept { return _data.mesh(); }

  };
}

#endif
