#include "msh/mesh_shape_interplay.hpp"
#include "apt/block.hpp"

namespace msh:: impl {
  // RATIONALE Due to the MIDWAY offset of dJ field, we define the NATIVE grid te be one shifted to the right by half spacing from the standard grid, therefore q_nat = q_std - 0.5. By design, dJ[i] is INSITU on the native grid.
  // The contributed cells in the native grid are [ INT_FLOOR( q_nat_min - sf.r) + 1,  INT_FLOOR( q_nat_max + sf.r) + 1 ). These are also the correct cell numbers in the standard grid
  // Now we have q0_std and q1_std, the final range of contributed cells is the union of individual ones

  template < int DGrid, typename Tq, typename ShapeF >
  constexpr auto deposit_range( const Tq& q0_std, const Tq& q1_std, const ShapeF& ) {
    static_assert( ShapeF::support() > 0 );
    apt::Index<DGrid> I_b; // index of the first contributed cell in the standard grid
    apt::Index<DGrid> I_e;

    // function of support = 1, jumps at support boundaries. The following assumes a left continuous convention, i.e. S(x) = lim S(x - \epsilon). In other words, nontrivial influence is only for left-open-right-closed ranges. It has to do with int(0) = 0.
    static_assert( ShapeF::support() > 1); // TODOL
    auto int_flr =
      []( auto q ) noexcept -> int {
        // Since min(q0_std) = 0.0 by design, min(q_nat) = -1 - r - 0.5, so q_nat + ( support + 3 ) / 2.0 >= 0. We will simply use (support + 1) as the shift
        return int(q + 1 + ShapeF::support() ) - ( 1 + ShapeF::support() );
      };

    for ( int i = 0; i < DGrid; ++i ) {
      auto a = q0_std[i];
      auto b = q1_std[i];

      // initially, a = q0_std, b = q1_std.
      b -= a;
      a -= 0.5;
      // now, a = q0_nat, b = dq_nat.
      a += ( b < 0 ? b : 0.0 );
      b = a + b * ( (b > 0.0) - ( b < 0.0) );
      // now a = min(q0_nat, q1_nat), b = max(q0_nat,q1_nat)

      I_b[i] = int_flr( a - ShapeF::support() / 2.0 ) + 1;
      I_e[i] = int_flr( b + ShapeF::support() / 2.0 ) + 1;
    }

    return apt::pair<decltype(I_b)>{ I_b, I_e };
  }

}

namespace msh {
  template < typename T >
  constexpr T Wesir( T sx0, T sx1, T sy0 = 1.0, T sy1 = 1.0 ) noexcept {
    return ( ( 2 * sx1 + sx0 ) * sy1 + ( sx1 + 2 * sx0 ) * sy0 ) / 6.0;
  }

  template < typename RealJ, int DField, int DGrid, typename ShapeF, typename U >
  void deposit ( field::Field<RealJ,DField,DGrid>& _J,
                 U charge_over_dt,
                 const ShapeF& shapef,
                 const apt::array<U,DField>& q0_std,
                 const apt::array<U,DField>& q1_std ) {
    static_assert( DField == 3 );
    static_assert( DGrid > 1 && DGrid < 4 );

    const auto[I_b, I_e] = impl::deposit_range<DGrid>( q0_std, q1_std, shapef );

    constexpr auto s =
      [] ( int i, auto q ) noexcept {
        return ShapeF()(i + 0.5 - q);
      };

    static_assert( DGrid == 2 ); // TODO DGrid == 3 not implemented yet

    // POLEDANCE use Block in the rest of code
    const auto& mesh = _J.mesh();
    if constexpr ( DGrid == 2 ) {
        for ( int n = 0; n < 2; ++n ) { // n is normal direction, 1 - n is transverse
          int n_stride = mesh.stride()[n];
          for ( int tr = I_b[1-n]; tr < I_e[1-n]; ++ tr ) {
            int trIdx = mesh.linear_index(apt::Longidx(1-n,tr));

            U Wform = Wesir( s(tr,q0_std[1-n]), s(tr,q1_std[1-n]) ) * charge_over_dt;
            U jsum = 0.0;
            for ( int i = I_e[n] - 1; i > I_b[n]; --i ) { // NOTE current at I_b[n] is zero by the normalization of shape function, so skip it
              jsum += s(i,q1_std[n]) - s(i,q0_std[n]);
              _J[n][i*n_stride + trIdx] += jsum * Wform;
            }
          }
        }

        for ( int j = I_b[1]; j < I_e[1]; ++j ) {
          for ( int i = I_b[0]; i < I_e[0]; ++i ) {
            _J[2]({i,j}) += (q1_std[2] - q0_std[2]) * Wesir( s(i,q0_std[0]), s(i,q1_std[0]), s(j,q0_std[1]), s(j,q1_std[1]) ) * charge_over_dt;
          }
        }

      }

  }

}
