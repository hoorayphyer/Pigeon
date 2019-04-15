#include "field/current_deposition.hpp"
#include "field/communication.hpp"
#include "apt/block.hpp"
#include "parallel/mpi++.hpp"

namespace field :: impl {
  // RATIONALE Due to the MIDWAY offset of dJ field, we define the NATIVE grid te be one shifted to the right by half spacing from the standard grid, therefore q_nat = q_std - 0.5. By design, dJ[i] is INSITU on the native grid.
  // The contributed cells in the native grid are [ INT_FLOOR( q_nat_min - sf.r) + 1,  INT_FLOOR( q_nat_max + sf.r) + 1 ). These are also the correct cell numbers in the standard grid
  // Now we have q0_std and q1_std, the final range of contributed cells is the union of individual ones

  template < int DGrid, typename Tq, typename ShapeF >
  constexpr auto deposit_range( const Tq& q0_std, const Tq& q1_std, const ShapeF& ) {
    static_assert( ShapeF::support() > 0 );
    apt::Index<DGrid> I_b; // index of the first contributed cell in the standard grid
    apt::Index<DGrid> extent;

    apt::foreach<0, DGrid>
      ( [] ( auto& i_b, auto& ext, auto a, auto b )
        noexcept {
          // initially, a = q0_std, b = q1_std.
          b -= a;
          a -= 0.5;
          // now, a = q0_nat, b = dq_nat.
          a += ( b < 0 ? b : 0.0 );
          b = a + b * ( (b > 0.0) - ( b < 0.0) );
          // now a = min(q0_nat, q1_nat), b = max(q0_nat,q1_nat)

          // function of support = 1, jumps at support boundaries. The following assumes a left continuous convention, i.e. S(x) = lim S(x - \epsilon). In other words, nontrivial influence is only for left-open-right-closed ranges. It has to do with int(0) = 0.
          static_assert( ShapeF::support() > 1); // TODOL

          auto int_flr =
            []( auto q ) noexcept {
              // Since min(q0_std) = 0.0 by design, min(q_nat) = -1 - r - 0.5, so q_nat + ( support + 3 ) / 2.0 >= 0. We will simply use (support + 1) as the shift
              constexpr auto shift = 1 + ShapeF::support();
              return int(q+shift) - shift; };

          i_b = int_flr( a - ShapeF::support() / 2.0 ) + 1 ;
          ext = int_flr( b + ShapeF::support() / 2.0 ) + 1 - i_b;
        },
        I_b, extent, q0_std, q1_std );

    return apt::pair<decltype(I_b)>{ I_b, extent };
  }

}

namespace field {
  template < typename T, int DField, int DGrid, typename ShapeF >
  Standard_dJ_Field<T,DField,DGrid,ShapeF>
  ::Standard_dJ_Field( apt::Index<DGrid> bulk_extent, const ShapeF& shapef )
    // NOTE minimum needed number of guards on one side is ( supp + 1 ) / 2 + 1
    : _data({ std::move(bulk_extent), ( ShapeF::support() + 3 ) / 2 }) {
    // enforce offset
    apt::array< offset_t, DGrid > offset{};
    apt::foreach<0,DGrid>
      ( []( auto& ofs ){ ofs = MIDWAY; }, offset );
    for ( int i = 0; i < DField; ++i )
      _data.set_offset( i, offset );
  }

  template < typename T >
  constexpr T Wesir( T sx0, T sx1, T sy0 = 1.0, T sy1 = 1.0 ) noexcept {
    return ( ( 2 * sx1 + sx0 ) * sy1 + ( sx1 + 2 * sx0 ) * sy0 ) / 6.0;
  }


  template < typename T, int DField, int DGrid, typename ShapeF >
  template < typename U >
  void Standard_dJ_Field<T,DField,DGrid,ShapeF>
  ::deposit ( U charge_over_dt, const apt::array<U,DField>& q0_std, const apt::array<U,DField>& q1_std ) {
    static_assert( DField == 3 );
    static_assert( DGrid > 1 && DGrid < 4 );

    constexpr auto shapef = ShapeF();
    const auto[I_b, extent] = impl::deposit_range<DGrid>( q0_std, q1_std, shapef );

    apt::array<U,3> W{};
    apt::array<U,DGrid> s0{};
    apt::array<U,DGrid> s1{};

    auto calc_s =
      [&shapef] ( auto& s, int i, auto q ) noexcept {
        s = shapef(i + 0.5 - q);
      };

    for ( auto I : apt::Block(extent) ) {
      I += I_b;
      // TODOL use normalization?? er ci xing zheng ze hua
      apt::foreach<0,DGrid> ( calc_s, s0, I, q0_std );
      apt::foreach<0,DGrid> ( calc_s, s1, I, q1_std );

      if constexpr ( DGrid == 2 ) {
          W[0] = ( s1[0] - s0[0] ) * Wesir( s0[1], s1[1] );
          W[1] = ( s1[1] - s0[1] ) * Wesir( s0[0], s1[0] );
          W[2] = (q1_std[2] - q0_std[2]) * Wesir( s0[0], s1[0], s0[1], s1[1] );
        } else {
        W[0] = ( s1[0] - s0[0] ) * Wesir( s0[1], s1[1], s0[2], s1[2] );
        W[1] = ( s1[1] - s0[1] ) * Wesir( s0[2], s1[2], s0[0], s1[0] );
        W[2] = ( s1[2] - s0[2] ) * Wesir( s0[0], s1[0], s0[1], s1[1] );
      }

      // NOTE: Before this line, there is no massive number of additions hence no loss of precision, which may only happen during the following +=. The fact that dj has higher precision takes care of all that.
      apt::foreach<0, DField> // NOTE it is DField here, not DGrid
        ( [&]( auto dj, auto w ) { dj(I) += w * charge_over_dt; } , _data, W ); // TODOL semantics on dj
    }

    return;
  }

}

namespace field {
  template < typename T, int DField, int DGrid, typename ShapeF >
  void Standard_dJ_Field<T,DField,DGrid,ShapeF>
  ::reduce( int chief, const mpi::Comm& intra ) {
    // TODO Opimize communication. Use persistent and buffer? NOTE one can use one chunk of memory for Jmesh so that only one pass of reduce is needed ?
    for ( int i = 0; i < DField; ++i )
      intra.reduce<mpi::IN_PLACE>( mpi::by::SUM, chief, _data[i].data().data(),  _data[i].data().size() );
  }
}

namespace field {
  template < typename T, int DField, int DGrid, typename ShapeF >
  Field<T,DField,DGrid>& Standard_dJ_Field<T,DField,DGrid,ShapeF>
  ::integrate( const mpi::CartComm& cart ) {
    auto& dJ = _data;

    const auto& mesh = dJ.mesh();
    for ( int i_dim = 0; i_dim < DGrid; ++i_dim ) { // NOTE it is DGrid not DField
      auto comp = dJ[i_dim]; // TODOL semantics
      for ( auto trI : mesh.project( i_dim, mesh.origin(), mesh.extent() ) ) {
        for ( int n = mesh.bulk_dim(i_dim) + mesh.guard() - 2; n > -mesh.guard() - 1; --n )
          comp[trI | n] += comp[trI | n+1];
      }
    }

    field::merge_guard_cells_into_bulk( dJ, cart );
    return dJ;
  }
}
