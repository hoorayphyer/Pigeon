#include "ParticleUpdater.h"
#include "AperParams.h"
#include "mpi_wrapper.h"
#include "Rng.h"

#include <chrono>
#include "Logger.h"
#include "InfoCollector.h"
#include "AperData.h"

#include "CurrentDepositer.cpp"
#include "ParticleInjector.cpp"
#include "ParticlePusher.cpp"
#include "ParticleSorter.h"
#include "ParticleCommunicator.h"
#include "PairProducer.h"
#include "AnnihManager.h"

#include "Domain.h"

using namespace std::chrono;

template< typename ...Rest >
inline double ShowDuration( high_resolution_clock::time_point& t_last,
                          const std::string& name_action,
                          const Rest&... rest ) {
  auto t_now = high_resolution_clock::now();
  auto dur = duration_cast<clock_cycle>( t_now - t_last );
  Logger::print(0, "------ Time for ", name_action, rest..., "is", dur.count(), clock_cycle_str,  ". ------");
  t_last = high_resolution_clock::now();
  return static_cast<double>( dur.count() );
}

template <CoordType Coord>
ParticleUpdater<Coord>::ParticleUpdater( const Dashboard& db, const AperParams& params, const MPICommunicator& comm, Domain* domain, std::unordered_map<ParticleType, ParticleSorter*>& pSorter, ParticleCommunicator*& ptc_comm )
  : _pSorter(pSorter), _ptc_comm(ptc_comm) {
  const auto& grid = params.grid;

  _sort2ndInterval = db.sort2ndInterval;
  _annih_pane = db.annih;

#if DEBUGLOG
  Logger::debug_print(debug_indent*2, "ParticlePusher");
#endif
  _pPusher = new ParticlePusher( grid, db.pusher );

#if DEBUGLOG
  Logger::debug_print(debug_indent*2, "CurrentDepositer");
#endif
  _pDepositer = new CurrentDepositer( grid, domain );

#if DEBUGLOG
  Logger::debug_print(debug_indent*2, "ParticleInjector");
#endif
  // initialize injectors at the specified boundary for all processes, regardless of whether the process is at the real boundary or not. This ensures code functionality when dynamic adjust is turned on
  for ( const auto& elm : db.boundaryCondition.ptcBC ) {
    auto bpos = elm.first;
    const auto& pbc = elm.second;
    if ( pbc.type == ParticleBCType::INJECT )
      _pInjector.emplace( bpos, new ParticleInjector(bpos, params.Q_e, pbc) );
  }
#if DEBUGLOG
  Logger::debug_print(debug_indent*2, "PairProducer");
#endif
  switch (db.pairProduce.pairScheme) {
  case PairProductionScheme::INSTANT :
    _pPairProducer = new PairProducer_Instant( db.pairProduce); break;
  case PairProductionScheme::PHOTON :
    _pPairProducer = new PairProducer_Photon( db.pairProduce); break;
  default :
    _pPairProducer = new PairProducer;
  }

#if DEBUGLOG
  Logger::debug_print(debug_indent*2, "AnnihManager");
#endif
  _pAnnih = new AnnihManager;
}

template <CoordType Coord>
ParticleUpdater<Coord>::~ParticleUpdater() {
  delete _pPusher;
  delete _pDepositer;
  for ( auto& ptr : _pInjector )
    delete ptr.second;
  delete _pPairProducer;
  delete _pAnnih;
}


// FIXME remove this. Only need ptcNumbers for inject
template < typename P >
void CountNumbers( const Particles<P>& particles, ScalarField<unsigned int>& ptcNum ) {
  ptcNum.assign(0);
  for ( unsigned int i = 0 ; i < particles.Number(); ++i ) {
    if ( particles.IsEmpty(i) ) continue;
    ptcNum.data() [ particles.PtcData()[i].cell ] += 1;
  }

}


template <CoordType Coord>
void ParticleUpdater<Coord>::Update(AperData& data, int timestep, const AperParams& params, const MPICommunicator& comm, Rng& rng, Domain* domain) {
  const auto dt = params.dt;
  const auto& grid = params.grid;
  const auto weightType = params.weightType;
  const auto& is_bdry = params.ens_specs.is_at_boundary;
  const auto& is_axis = params.ens_specs.is_axis;

  auto& info = InfoCollector::Instance().contents;
  auto& ssProxy = InfoCollector::Instance().ssProxy;

  comm.ensemble().barrier(); // for better timing on EnsBroadCastFields
#if DEBUGLOG
  Logger::debug_print(debug_indent*3, "EnsBroadCastFields");
#endif
  auto t0 = high_resolution_clock::now();
  auto t_last = t0;
  domain -> EnsBroadcastFields( data.Efield, data.Bfield );
  info.t_ens_bcastEB = ShowDuration( t_last, "Ens broadcast E&B" );

  //------ Update Particle x and p ------
  // NOTE by this order, x is half a time step ahead of p. This makes it simple to export J by species because at that time x is not updated by dx yet.
  // info.t_push_qcarriers = 0.0;
  for ( auto& elm : data.particles ) {
    auto ptcType = elm.first;
    auto& ptcs = elm.second;
    // FIXME TODO should deal with rebase momentum in pusher
#if DEBUGLOG
    Logger::debug_print(debug_indent*3, "Push", PtcType2Str(ptcType));
#endif
    if ( is_charged(ptcType) )
      _pPusher -> Push<Coord>( QC(ptcs), dt, data.Efield, data.Bfield, params );
    else
      _pPusher -> Push<Coord>( QN(ptcs), dt, params );
    // info.t_push_qcarriers += ShowDuration( t_last, "push" );
  }

  //------ Sort and SendRecv ------
  // info.t_sort1_qcs = 0.0;
  // info.t_sendrecv_qcs = 0.0;
  // info.t_sort2_ptc = 0.0;
  for ( auto& elm : data.particles ) {
    auto ptcType = elm.first;
    auto& ptcs = elm.second;
    auto& sorter = *_pSorter.at(ptcType);
    auto ptcStr = PtcType2Str(ptcType);

#if DEBUGLOG
    Logger::debug_print(debug_indent*3, "Sort1", ptcStr );
#endif
    std::visit( [&]( const auto& x )
                { sorter.Sort( x, x.Number(), grid, true); },
      ptcs );
    // info.t_sort1_qcs += ShowDuration( t_last, "sorting", ptcs.NameStr(),"1st time" );

#if DEBUGLOG
    Logger::debug_print(debug_indent*3, "SendRecvParticles", ptcStr );
#endif
    auto& ptc_comm = *_ptc_comm;
    std::visit( [&]( const auto& x )
                { ptc_comm.SendRecvParticles( timestep, x, sorter, comm, params ); },
      ptcs );
    // info.t_sendrecv_qcs += ShowDuration( t_last, "sendrecv", ptcs.NameStr() );
    if ( timestep % _sort2ndInterval == 0 ) {
#if DEBUGLOG
      Logger::debug_print(debug_indent*3, "Sort2", ptcStr);
#endif
      std::visit( [&]( const auto& x )
                  { sorter.Sort( ptcs, ptcs.Number(), grid, false);},
                  ptcs );
      // info.t_sort2_ptc += ShowDuration( t_last, "sorting", particles.NameStr(),"2nd time" );
    }

  }

  info.num_annihPairs = 0;
  const bool is_annih_ep = _annih_pane.matter.isOn && (timestep % _annih_pane.matter.interval == 0);
  // FIXME TODO use reduce-scatter for annih_nums
  // FIXME TODO limit outstanding reqs in EnsReduceFields in Domian.
  if ( is_annih_ep ) {
#if DEBUGLOG
    Logger::debug_print(debug_indent*3, "AnnihMark ep");
#endif
    _pAnnih -> Mark( data.get<ParticleType::ELECTRON>(), data.get<ParticleType::POSITRON>(), dt, grid, comm.ensemble(), _annih_pane.matter );
    info.t_AnnihMark_ep = ShowDuration( t_last, "AnnihMark ep" );
  }

#if DEBUGLOG
  Logger::debug_print(debug_indent*3, "DepositCurrent");
#endif
  info.t_deposit_mpi = 0.0; // DepositCurrent will add to this variable
  // FIXME ptc interface
  _pDepositer->DepositCurrent<Coord>( data, dt, weightType, comm, is_bdry, is_axis);
  info.t_deposit = ShowDuration( t_last, "depositing current" );


  // ------ birth and death of particles ------
#if DEBUGLOG
  Logger::debug_print(debug_indent*3, "Annihilate particles");
#endif
  if ( is_annih_ep )
    _pAnnih -> Annihilate( QC(data.particles.at(ParticleType::ELECTRON)), QC(data.particles.at(ParticleType::POSITRON)) );
  _pAnnih -> Annihilate( QN( data.particles.at(ParticleType::PHOTON) ), grid, _annih_pane.photon );

#if DEBUGLOG
  Logger::debug_print(debug_indent*3, "InjectPairs");
#endif
  info.num_injectedPairs = 0;
  for ( auto& inj : _pInjector ) {
    if ( params.ens_specs.is_at_boundary[inj.first] ) // only inject at the real boundary
      // FIXME ptc interface
      inj.second -> InjectPairs<Coord>( timestep, dt, data, weightType, comm, rng );
  }
  info.t_injection = ShowDuration( t_last, "injection" );

#if DEBUGLOG
  Logger::debug_print(debug_indent*3, "ProducePairsAndPhotons");
#endif
  // FIXME ptc interface
  _pPairProducer->ProducePairsAndPhotons( timestep, data, dt, grid, rng );
  info.t_pairPhotonCreation = ShowDuration( t_last, "pair production" );


  //------ Wrap up ----------
  // FIXME
  // auto f_count_tracked =
  //   [] ( const auto& particles ) {
  //     int sum = 0;
  //     for ( int i = 0; i < particles.Number(); ++i ) {
  //       if ( !particles.IsEmpty(i) && ParticleTracker::IsTracked( particles.PtcData()[i].track_id ) )
  //         sum++;
  //     }
  //     return sum;
  //   };
  // for ( const auto& elm : data.particles ) {
  //   auto ptcType = elm.first;
  //   const auto& particles = elm.second;
  //   info.num_tot[static_cast<int>(ptcType)] = particles.Number();
  //   info.num_tracked[static_cast<int>(ptcType)] = f_count_tracked(particles);
  // }

  // info.num_tot[static_cast<int>(ParticleType::PHOTON)] = data.photons.Number();
  // info.num_tracked[static_cast<int>(ParticleType::PHOTON)] = f_count_tracked(data.photons);
  // info.size_pcEvents = data.pairCreationTracker.GetData().size();

  // for ( auto it = data.particles.cbegin(); it != data.particles.cend(); ++it ) {
  //   auto ptcType = it->first;
  //   CountNumbers( data.particles.at(ptcType), data.ptcNumbers.at(ptcType) );
  // }

  auto dur = duration_cast<clock_cycle>(high_resolution_clock::now() - t0);
  Logger::print(0, "------ Time for a particle update step in total is", dur.count(), clock_cycle_str, ". ------");
  info.t_updateParticle = dur.count();

}
