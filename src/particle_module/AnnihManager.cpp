#include "AnnihManager.h"
#include "mpi_wrapper.h"
#include "InfoCollector.h"

auto GetAnnihNums( Scalar dt, const ScalarField<unsigned int>& num_e, const ScalarField<unsigned int>& num_p, const MPIEnsembleCommunicator& ensemble, const DBPane_Annihilation::Matter& pane ) {
  const auto& grid = num_e.grid();
  ScalarField<unsigned int> annih_nums (grid); // OK to not zero out

  if ( ensemble.is_root() ) {
    ScalarField<unsigned int> num_tot_e(grid), num_tot_p(grid);
    ensemble.reduce( num_e.ptr(), num_tot_e.ptr(), grid.size(), MPI_SUM, ensemble.root() );
    ensemble.reduce( num_p.ptr(), num_tot_p.ptr(), grid.size(), MPI_SUM, ensemble.root() );

    auto& count = InfoCollector::Instance().contents.num_annihPairs;
    count = 0;
    Scalar annih_percent = std::min(1.0, pane.rate * pane.interval * dt ); // NOTE the pane.interval
    for ( int cell = 0; cell < grid.size(); ++cell ) {
      auto pos = grid.pos_3d( cell, {0,0,0} );
      if ( !pane.annihilable( pos.x, pos.y, pos.z ) ) continue;

      int num_pair = std::min( num_tot_e.data()[cell], num_tot_p.data()[cell] );
      if ( num_pair > pane.num_pairs_target ) {
        annih_nums.data()[cell] = num_pair * annih_percent;
        count += annih_nums.data()[cell];
      }
    }
  } else {
    ensemble.reduce( num_e.ptr(), nullptr, grid.size(), MPI_SUM, ensemble.root() );
    ensemble.reduce( num_p.ptr(), nullptr, grid.size(), MPI_SUM, ensemble.root() );
  }

  return annih_nums;
}

void MarkImpl( Particles<QCarrier>& ptcs, ScalarField<unsigned int>& annih_nums ) {
  for ( int i = 0; i < ptcs.Number(); ++i ) {
    if ( ptcs.IsEmpty(i) ) continue;

    auto& ptc = ptcs.PtcData()[i];
    if ( 0 == annih_nums.data()[ptc.cell] ) continue;

    set_bit( ptc.flag, ParticleFlag::annihilate );

    // move particles to the center of the cell. It is OK to do it in all 3 dimensions regardless of the real dimension of the simulation
    ptc.dx.x = 0.5 - ptc.x.x;
    ptc.dx.y = 0.5 - ptc.x.y;
    ptc.dx.z = 0.5 - ptc.x.z;

    --annih_nums.data()[ptc.cell];
  }
}

// FIXME duplicates
template < typename P >
void CountNumbers( const Particles<P>& particles, ScalarField<unsigned int>& ptcNum ) {
  ptcNum.assign(0);
  for ( unsigned int i = 0 ; i < particles.Number(); ++i ) {
    if ( particles.IsEmpty(i) ) continue;
    ptcNum.data() [ particles.PtcData()[i].cell ] += 1;
  }

}

void AnnihManager::MarkPairs( Particles<QCarrier>& electrons, Particles<QCarrier>& positrons, Scalar dt, const Grid& grid, const MPIEnsembleCommunicator& ensemble, const DBPane_Annihilation::Matter& pane ) {
  ScalarField<unsigned int> num_e(grid), num_p(grid);
  CountNumbers(electrons, num_e);
  CountNumbers(positrons, num_p);

  // annih_nums on root is regarded as a const. annih_nums on nonroot is not significant
  auto annih_nums = GetAnnihNums ( dt, num_e, num_p, ensemble, pane );
  const int rank = ensemble.rank();

  auto* ep[2] = { &electrons, &positrons };
  for ( auto& p : ep ) {
    if ( ensemble.is_root() ) {
      ScalarField<unsigned int> nums( annih_nums );
      MarkImpl( *p, nums );
      if ( rank != ensemble.size() - 1 )
        ensemble.send( rank + 1, 147, nums.ptr(), nums.gridSize() );
    } else {
      ensemble.recv( rank-1, 147, annih_nums.ptr(), annih_nums.gridSize() );
      MarkImpl( *p, annih_nums );
      if ( rank != ensemble.size() - 1 )
        ensemble.send( rank + 1, 147, annih_nums.ptr(), annih_nums.gridSize() );
    }
  }

}

void AnnihManager::Annihilate( Particles<QCarrier>& electrons, Particles<QCarrier>& positrons ) {
  auto f = [] ( auto& ptcs ) {
             for ( unsigned int i = 0; i < ptcs.Number(); ++i ) {
               if ( check_bit( ptcs.PtcData()[i].flag, ParticleFlag::annihilate ) )
                 ptcs.Erase( i, 1 );
             }
           };
  f(electrons);
  f(positrons);
}

void AnnihManager::Annihilate( Particles<QNeutral>& photons, const Grid& grid, const DBPane_Annihilation::Photon& pane ) {
  for ( unsigned int i = 0; i < photons.Number(); ++i ) {
    if ( photons.IsEmpty(i) ) continue;

    auto& ptc = photons.PtcData()[i];
    auto pos = grid.pos_particle( ptc.cell, ptc.x );
    if ( pane.annihilable( pos.x, pos.y, pos.z ) )
      photons.Erase( i, 1 );
  }
}
