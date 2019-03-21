#include "particle/updater.hpp"

#include "particle/pusher.hpp"
#include "particle/pair_producer.hpp"
#include "particle/properties.hpp"

#include "field/field_shape_interplay.hpp"
#include <algorithm>

#include "dynamic_variables.hpp"
#include "parameters.hpp"
#include "kernel/coordinate.hpp"

#include "kernel/grid.hpp"
#include "kernel/shapef.hpp"

namespace particle::impl {
  template < species sp, PairScheme pair_scheme, knl::coordsys CS,
             typename DynaVars_t, typename JmeshField, typename Migrators, typename Params_t, typename ShapeF, typename Rng, typename Borders >
  void update_species( DynaVars_t& dvars, JmeshField& Jmesh, Migrators& migrators, const Params_t& params, const ShapeF& shapef, Rng& rng, const Borders& borders ) {
    if ( dvars[sp].size() == 0 ) return; // TODO this is wrong, unused particle != having zero particles

    auto dt = params.dt;

    for ( auto& ptc : dvars[sp] ) { // TODO check ptc is proxy
      if( ptc.is(flag::empty) ) continue;

      if constexpr ( is_charged<sp> ) {
          auto E = field::interpolate(dvars.E, ptc.q, shapef );
          auto B = field::interpolate(dvars.E, ptc.q, shapef );
          auto&& dp = update_p<sp>( ptc, dt, E, B );

          // TODO make sure the newly added particles will not participate in the loop
          if constexpr ( pair_scheme != PairScheme::Disabled && is_radiative<sp> ) {
            auto Rc = calc_Rc( dt, ptc.p, std::move(dp) );
            auto gamma = std::sqrt( is_massive<sp> + apt::sqabs(ptc.p) );
            if ( !is_productive_lepton( ptc, gamma, Rc, rng ) )
              goto end_of_pair_produce;

            if constexpr ( pair_scheme == PairScheme::Photon ) {
              produce_photons( std::back_inserter( dvars.photons ),
                               ptc, gamma, Rc );
            } else if ( pair_scheme == PairScheme::Instant ) {
              instant_produce_pairs( std::back_inserter( dvars.electrons ),
                                     std::back_inserter( dvars.positrons ),
                                     ptc, gamma, Rc );
            }
          }
          end_of_pair_produce:;

        }
      else if ( sp == species::photon && pair_scheme == PairScheme::Photon ) {
        photon_produce_pairs( std::back_inserter( dvars.electrons ),
                              std::back_inserter( dvars.positrons ),
                              ptc );
      }

      // NOTE q is updated, starting from here, particles may be in the guard cells.
      auto&& dq = update_q<sp,CS>( ptc, dt );
      // TODO pusher handle boundary condition
      if constexpr ( is_charged<sp> )
                     // NOTE: we use -dq here, so the deposited is TdJ ( time reversal of dJ ). This reversal will be cancelled by a reversed integration in integrate_TdJ due to our choice of indexing.
        field::deposit_dJ( Jmesh, charge_x<sp>, ptc.q(), -std::move(dq), shapef );

      if ( is_migrate( ptc.q(), borders ) ) {
        migrators.emplace_back(std::move(ptc));
      }

    }

    // TODO sort particle array periodically
  }

  template < PairScheme PS, knl::coordsys CS, int I = 0, typename... Args >
  inline void iterate_species( Args&&... args ) {
    if constexpr ( I == 4 ) return; // TODOL loop over all species. Would void_t help?
    else {
      // update_species<static_cast<species>(I), PS, CS>( std::forward<Args>(args)... );
      return iterate_species<PS, CS, I+1>( std::forward<Args>(args)... );
    }
  }

}

namespace particle::impl {
  template < typename TdJ_t >
  void integrate_TdJ( TdJ_t& TdJ, const mpi::CartComm& cart ) {
    // NOTE
    // - we assum TdJ is offset at MIDWAY in all directions
    // - the boundary conditions are J =0 at both ends in a direcction
    // - reminder: TdJ here is the time reversal of dJ
    // - reminder: indexing is such that 0 corresponds to the first cell in bulk
    // - reminder: guard on dJ is shape::support / 2 + 1
    // NOTE By convention of our indexing, J[i+1] = J[i] + dJ[i] = J[i] - TdJ[i], which is not in-place doable. So we use J[i] = J[i+1] + TdJ[i], i.e. we integrate from upper bound backwards. Given the boundary conditions, J[bulk_dim - guard] can be found locally. From there, one can find all J_i down through J[guard]. NOTE J[guard] can be found two ways, they should be consistent thanks to charge conservation of Esirkepov algorithm. See below for how to leverage this on achieving complete in-place integration.
    // NOTE: After merging cells, TdJ values in those untreated cells are also ready. In the lower end, one can continue J[i] = J[i+1] + TdJ[i] in place. However, at the upper end, one will need J[i+1] = J[i] - TdJ[i], which breaks inplace-ness. More importantly, the TdJ value of the very first cell ( call this cell FB ) counted from the back is overwritten with the J[FB], but the original TdJ[FB] is needed again to find J[FB + 1] = J[FB] - TdJ[FB]. Luckily charge conservation makes some information redundant. One can find that TdJ[-1] is not needed during merging cells ( because J[0] is found from J[1] + TdJ[0], while J[-1], which resides on the neighboring process, is found from J[-2] - TdJ[-2] ). So we will use it to store J[FB] at the appropriate time.

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

  }
}

namespace particle {
  template < typename Real, int DGrid, int DPtc, typename state_t,
             knl::shape Shape, PairScheme pair_scheme, knl::coordsys CS
             >
  void Updater< Real, DGrid, DPtc, state_t, Shape, pair_scheme, CS >
  ::operator() ( int timestep, DynaVars_t& dvars,
                 const Params_t& params,
                 const std::optional<mpi::CartComm>& cart,
                 const Ensemble<DGrid>& ensemble ) {

    // TODO borders at real physical boundary need to be specified still, because now mesh controls margin
    apt::array< apt::pair<Real>, DGrid > borders; // TODO borders not done yet

    apt::foreach<0, 3>
      ( []( auto comp ) { // returns a proxy
          for ( auto& elm : comp.data() ) elm = 0.0;
        }, _Jmesh );

    // TODOL reduce number of communications?
    for ( int i = 0; i < 3; ++i )
      ensemble.intra.broadcast( ensemble.chief, dvars.E[i].data().data(), dvars.E[i].data().size() );
    for ( int i = 0; i < 3; ++i )
      ensemble.intra.broadcast( ensemble.chief, dvars.B[i].data().data(), dvars.B[i].data().size() );

    impl::iterate_species<pair_scheme, CS>( dvars, _Jmesh, _migrators, params,
                                            knl::shapef_t<Shape>(), _rng, borders );

    migrate( _migrators, _intercomms, borders, timestep );
    for ( auto&& ptc : _migrators ) {
      dvars.particles[ptc.template get<species>()].push_back( std::move(ptc) ); // TODO check this
    }
    _migrators.resize(0);

    // TODO do this before J is actually used
    {
      const auto& grid = dvars.E.mesh().bulk();

      apt::foreach<0, DGrid>
        ( [&]( auto comp, const auto& g ) {
            auto tmp = -1 * params.e * g.delta() / params.dt;
            for ( auto& elm : comp.data() ) elm *= tmp;
          }, _Jmesh, grid );

      if constexpr ( DGrid == 2 ) {
        auto tmp = params.e / params.dt;
        for ( auto& elm : _Jmesh[2].data() ) elm *= tmp;
      }

      // TODO NOTE one can use one chunk of memory for Jmesh so that only one pass of reduce is needed
      for ( int i = 0; i < 3; ++i )
        ensemble.intra.reduce<by::SUM, mpi::IN_PLACE>( ensemble.chief, _Jmesh[i].data().data(),  _Jmesh[i].data().size() );

      if ( cart ) {
        impl::integrate_TdJ( _Jmesh, *cart );
      }

    }
  }

}

#include "traits.hpp"
using namespace traits;
namespace particle {
  template class Updater< real_t, DGrid, DPtc, ptc_state_t, shape,
                          pair_produce_scheme, coordinate_system>;
}
