#include "field/field_shape_interplay.hpp"
#include "field/communication.hpp"
#include "apt/block.hpp"
#include "apt/vec.hpp"
#include "apt/type_traits.hpp"
#include "parallel/mpi++.hpp"
#include <tuple>

namespace esirkepov {
  template < typename T >
  inline T Wesir_2D( T sx0, T sx1, T sy0, T sy1 ) noexcept {
    return ( ( 2 * sx1 + sx0 ) * sy1 + ( sx1 + 2 * sx0 ) * sy0 ) / 6.0;
  }

  template < typename T >
  inline T Wesir( T sx0, T sx1, T sy0, T sy1, T sz0, T sz1 ) noexcept {
    return (sx1 - sx0) * Wesir_2D(sy0, sy1, sz0, sz1);
  }

  template < typename T, int DGrid, typename ShapeF >
  constexpr apt::array<T,3> calcW( const apt::Index<DGrid>& index,
                                   const apt::array<T, DGrid>& sep0_b,
                                   const apt::array<T, DGrid>& sep1_b,
                                   const ShapeF& shapef ) noexcept {
    T W[3]{};

    // TODOL use normalization?? er ci xing zheng ze hua
    T s0[DGrid]{};
    T s1[DGrid]{};

    auto f_calc_s0s1
      = [&shapef] ( auto& s0_i, auto& s1_i, auto i, auto sep0_i, auto sep1_i ) {
          s0_i = shapef( i + sep0_i );
          s1_i = shapef( i + sep1_i );
        };

    apt::foreach<0,DGrid>( f_calc_s0s1, s0, s1, index, sep0_b, sep1_b );

    if constexpr ( DGrid == 2 ) {
        W[0] = Wesir( s0[0], s1[0], s0[1], s1[1], 1.0, 1.0 );
        W[1] = Wesir( s0[1], s1[1], 1.0, 1.0, s0[0], s1[0] );
        W[2] = Wesir( 0.0, 1.0, s0[0], s1[0], s0[1], s1[1] );
      } else {
      W[0] = Wesir( s0[0], s1[0], s0[1], s1[1], s0[2], s1[2] );
      W[1] = Wesir( s0[1], s1[1], s0[2], s1[2], s0[0], s1[0] );
      W[2] = Wesir( s0[2], s1[2], s0[0], s1[0], s0[1], s1[1] );
    }

    return {{ W[0], W[1], W[2] }};
  }

}

namespace field :: impl {

  // RATIONALE: assume the offset = 0.5 in the dimension. The native grid is one offset by 0.5 with respect to the original one.
  // - q is relative position in the native grid. q starts at 0.0
  // - the contributing cells in the native grid are [ int(q - sf.r) + 1,  int(q + sf.r) + 1 )
  // - index_original = index_native + ( q - int(q) >= 0.5 )
  // Now we have q0 and q1, the final range of contributing cells is the union of individual ones
  template < typename Vec_q0_abs, typename Vec_q1_abs, typename Grid, typename ShapeF >
  constexpr auto set_up_for_depositing_dJ( const Vec_q0_abs& q0_abs, const Vec_q1_abs& q1_abs, const Grid& grid, const ShapeF& ) {
    constexpr int DGrid = apt::ndim_v<Grid>;
    using T = apt::most_precise_t<apt::element_t<Vec_q0_abs>, apt::element_t<Vec_q1_abs>>;
    apt::Index<DGrid> I_b;
    apt::Index<DGrid> extent;
    apt::array<T, DGrid> sep0_b; // separation 0 is defined to be i_cell_beginning - q0_ptc
    apt::array<T, DGrid> sep1_b;

    apt::foreach<0, DGrid>
      ( [] ( auto& i_b, auto& ext, auto& sp0, auto& sp1,
             auto qmin, auto qmax, const auto& g1d )
        noexcept {
          // initially, qmin = q0_abs, qmax = q1_abs. Find dq_abs = q1_abs - q0_abs and store it in qmax
          qmax -= qmin;
          qmax /= g1d.delta();
          qmin = ( qmin - g1d.lower() ) / g1d.delta() - 0.5;
          // now, qmin = q0_native, qmax = dq_native.
          sp0 = - qmin;
          sp1 = sp0 - qmax;
          qmin += ( ( qmax < 0 ? qmax : 0.0 ) - ShapeF::support / 2.0 );
          qmax = qmin + qmax * ( (qmax > 0.0) - ( qmax < 0.0) ) + ShapeF::support;
          i_b = int( qmin ) + 1 ;
          ext = int( qmax ) + 1 - i_b;
          sp0 += i_b;
          sp1 += i_b;
        },
        I_b, extent, sep0_b, sep1_b,
        q0_abs, q1_abs, grid );

    return std::make_tuple( I_b, extent, sep0_b, sep1_b );
  }

}

namespace field {
  template < typename T, int DField, int DGrid >
  dJ_Field<T,DField,DGrid>::dJ_Field( const Mesh<T,DGrid>& mesh )
    : _data(mesh) {
    // TODO need shape support + 1
    // enforce offset
    apt::array< offset_t, DGrid > offset{};
    apt::foreach<0,DGrid>
      ( []( auto& ofs ){ ofs = MIDWAY; }, offset );
    for ( int i = 0; i < DField; ++i )
      _data.set_offset( i, offset );
  }

  template < typename T, int DField, int DGrid >
  template < typename Vec_q0, typename Vec_q1, typename ShapeF >
  void dJ_Field<T,DField,DGrid>::deposit ( T charge_over_dt,
                                           const apt::VecExpression<Vec_q0>& q0_abs,
                                           const apt::VecExpression<Vec_q1>& q1_abs,
                                           const ShapeF& shapef ) {
    static_assert( DField == 3 );
    static_assert( DGrid > 1 && DGrid < 4 );

    // NOTE q0 and q1 are switched positions so that the deposited is the time reversal of dJ, or TdJ. The reason for this is elaborated in integrate.
    // NOTE: we use -dq here, so the deposited is TdJ ( time reversal of dJ ). This reversal will be cancelled by a reversed integration in integrate_TdJ due to our choice of indexing.
    const auto[I_b, extent, sep0_b, sep1_b] = set_up_for_depositing_dJ( q1_abs, q0_abs, _data.mesh().bulk(), shapef );

    for ( const auto& I : apt::Block(extent) ) {
      auto W = esirkepov::calcW( I, sep0_b, sep1_b, shapef );
      if constexpr ( DGrid == 2 ) W[2] *= (q1_abs[2] - q0_abs[2]); // NOTE here the forward time dq is used because there is no integration on dJ[2] in 2D case. TODO check
      apt::foreach<0, DField> // NOTE it is DField here, not DGrid
        ( [&]( auto& dj, auto w ) { dj(I_b + I) += w * charge_over_dt; } , _data, W );
    }

    // NOTE in eq.36 in Esirkepov, V_z is really del_z / dt, where del_z should
    // be interpretted as the COORDINATE displacement of the particle and dt is
    // the timestep. It is wrong to use the physical velcocity v_ptc_z here,
    // which is what we did before in LogSpherical where there was always a
    // confusion about where to base the v_ptc_z.

    return;
  }

}

namespace field {
  template < typename T, int DField, int DGrid >
  void dJ_Field<T,DField,DGrid>::reduce( int chief, const mpi::Comm& intra ) {
    // TODO NOTE one can use one chunk of memory for Jmesh so that only one pass of reduce is needed ?
    for ( int i = 0; i < DField; ++i )
      intra.reduce<mpi::by::SUM, mpi::IN_PLACE>( chief, _data[i].data().data(),  _data[i].data().size() );
  }
}

namespace field {
  template < typename T, int DField, int DGrid >
  const Field<T,DField,DGrid>& dJ_Field<T,DField,DGrid>::integrate( const mpi::CartComm& cart ) {
    // NOTE
    // - we assum TdJ is offset at MIDWAY in all directions
    // - the boundary conditions are J =0 at both ends in a direcction
    // - reminder: TdJ here is the time reversal of dJ
    // - reminder: indexing is such that 0 corresponds to the first cell in bulk
    // - reminder: guard on dJ is shape::support / 2 + 1
    // NOTE By convention of our indexing, J[i+1] = J[i] + dJ[i] = J[i] - TdJ[i], which is not in-place doable. So we use J[i] = J[i+1] + TdJ[i], i.e. we integrate from upper bound backwards. Given the boundary conditions, J[bulk_dim - guard] can be found locally. From there, one can find all J_i down through J[guard]. NOTE J[guard] can be found two ways, they should be consistent thanks to charge conservation of Esirkepov algorithm. See below for how to leverage this on achieving complete in-place integration.
    // NOTE: After merging cells, TdJ values in those untreated cells are also ready. In the lower end, one can continue J[i] = J[i+1] + TdJ[i] in place. However, at the upper end, one will need J[i+1] = J[i] - TdJ[i], which breaks inplace-ness. More importantly, the TdJ value of the very first cell ( call this cell FB ) counted from the back is overwritten with the J[FB], but the original TdJ[FB] is needed again to find J[FB + 1] = J[FB] - TdJ[FB]. Luckily charge conservation makes some information redundant. One can find that TdJ[-1] is not needed during merging cells ( because J[0] is found from J[1] + TdJ[0], while J[-1], which resides on the neighboring process, is found from J[-2] - TdJ[-2] ). So we will use it to store J[FB] at the appropriate time.
    auto& TdJ = _data;

    const auto& mesh = TdJ.mesh();
    for ( int i_dim = 0; i_dim < 3; ++i_dim ) {
      apt::Index<DGrid> I_b = mesh.origin();
      auto ext_trans = mesh.extent();
      ext_trans[i_dim] = 1;

      auto& comp = TdJ[i_dim]; // TODOL semantics

      // get range of cells whose J can be found locally. This covers ( J[guard-1], J[bulk_dim - guard] ].
      int iback_b = mesh.bulk()[i_dim].dim() - mesh.guard();

      // locally scan
      for ( auto I : apt::Block(ext_trans) ) {
        I += I_b;
        // first, store TdJ[bulk_dim - guard] and find J[bulk_dim - guard]. Set TdJ[-1] = 0 to avoid corrupting the stored value during merging guard cells
        std::swap( comp(iback_b), comp( iback_b + mesh.guard() - 1 ) );
        std::swap( comp(-1), comp( -mesh.guard() ) );
        comp( - mesh.guard() ) = 0.0; // TODO fix this

        for ( int i = iback_b + 1; i < iback_b + 2 * mesh().guard; ++i )
          comp( iback_b ) += comp( i );

        // then, perform scan
        for ( int i = iback_b - 1; i > mesh.guard() - 1; --i ) // NOTE --i
          comp( i ) += comp( i+1 );
      }
    }

    field::merge_guard_cells_into_bulk( TdJ, cart );

    // boundaries
    for ( int i_dim = 0; i_dim < 3; ++i_dim ) {
      apt::Index<DGrid> I_b{};
      auto ext_trans = mesh.bulk().extent();
      ext_trans[i_dim] = 1;

      auto& comp = TdJ[i_dim]; // TODOL semantics
      int iback = mesh.bulk()[i_dim].dim() - 1;

      for ( auto I : apt::Block(ext_trans) ) {
        // finish lower end
        for ( int i = mesh.guar() - 1; i > -1; --i ) // NOTE --i
          comp( i ) += comp( i+1 );

        // finish upper end
        // first find J[bulk_dim - 1]
        for ( int i = iback - mesh.guard() + 1; i < iback; ++i )
          comp( iback ) += comp( i );

        for ( int i = iback - 1; i > iback - mesh.guard() + 1; --i ) // NOTE ++i
          comp( i ) += comp( i+1 );
      }
    }

    apt::foreach<0, DGrid> // NOTE it's DGrid not DField
      ( [&]( auto comp, const auto& g ) { // comp is proxy
          auto tmp = -g.delta();
          for ( auto& elm : comp.data() ) elm *= tmp;
        }, TdJ, mesh.bulk() );

    return TdJ;
  }
}
