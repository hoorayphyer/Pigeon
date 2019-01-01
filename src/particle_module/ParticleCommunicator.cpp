#include "ParticleCommunicator.h"
#include "AperParams.h"
#include "mpi_wrapper.h"
#include "ParticleSorter.h"
#include "InfoCollector.h"
#include "EnsParticleGatherer.h"

// TagToTile and TileToTag functions depend strongly on how this enum is defined.
enum CommTagsParticles {
    FRONTBACK_DOWN_LEFT=0, FRONTBACK_DOWN, FRONTBACK_DOWN_RIGHT,
    FRONTBACK_LEFT, FRONTBACK, FRONTBACK_RIGHT,
    FRONTBACK_UP_LEFT, FRONTBACK_UP, FRONTBACK_UP_RIGHT,

    DOWN_LEFT, DOWN, DOWN_RIGHT,
    LEFT, CENTER, RIGHT,
    UP_LEFT, UP, UP_RIGHT,

    INVALID=127
};

inline int BufferVecToBufferIdx ( const Vec3<int>& bufferVec ) {
  // front and back buffers are degenerate
  int frontback = bufferVec.z == 0 ? 0 : -1;
  return bufferVec.x + bufferVec.y * 3 + frontback * 9 + 13;
}

inline int BufferVecToTile( const Vec3<int>& bufferVec) {
  return bufferVec.x + bufferVec.y*3 + bufferVec.z*9 + 13;
}

ParticleCommunicator::ParticleCommunicator(const MPICommunicator& comm) {
  // Initialize particle buffers
  // FIXME: How to find a more appropriate number?
  _sizeParticleBuffer = 100000;

  for(int i=0; i<COUNT_COMMTAGSPARTICLES; i++){
    _bufferParticles[i].Resize(_sizeParticleBuffer);
    _bufferParticles[i].Initialize();
    _bufferPhotons[i].Resize(_sizeParticleBuffer);
    _bufferPhotons[i].Initialize();
  }

  _ens_gatherer = new EnsParticleGatherer( comm.is_primary(), comm.ensemble().size() );
}

ParticleCommunicator::~ParticleCommunicator() {
  delete _ens_gatherer;
}

template <>
inline std::array<Particles<QCarrier>,COUNT_COMMTAGSPARTICLES>&
ParticleCommunicator::GetSendBuffer<QCarrier> () { return _bufferParticles; }

template <>
inline std::array<Particles<QNeutral>,COUNT_COMMTAGSPARTICLES>&
ParticleCommunicator::GetSendBuffer<QNeutral> () { return _bufferPhotons; }

template< typename P >
void ParticleCommunicator::SendRecvDirectional(Particles<P> &particles, const ParticleSorter& sorter, const MPICommunicator& comm, const Vec3<int>& sendVec, int direction, int neighbor_left, int neighbor_right){
    // check sendVec validity
    assert(sendVec[direction]!=0);
    for( int i=2; i>direction; i-- ) {
        assert(sendVec[i]==0);
    }

    bool is_send_right = sendVec[direction] > 0;
    int tile = BufferVecToTile( sendVec );
    int dest_rank = is_send_right? neighbor_right : neighbor_left;
    int src_rank = is_send_right ? neighbor_left : neighbor_right;
    const std::vector<int>& numInTile = sorter.getNumInTile();
    const std::vector<int>& scanTile = sorter.getScanTile();

    //prepare send and recv buffers
    auto& sendbuffer = GetSendBuffer<P>()[ BufferVecToBufferIdx(sendVec) ];
    assert( sendbuffer.Number() == 0 );

    if ( dest_rank != NEIGHBOR_NULL ) {
     sendbuffer.CopyFrom( particles, numInTile[tile], scanTile[tile], sendbuffer.Number() );

     if ( comm.ensemble().size() > 1 ) {
       auto& g_ensRecvBuf = _ens_gatherer->EnsBuffer<P>();
       // add particles gathered from replicas into sendbuffer as well
       int ens_buffer_idx = sendVec.x + sendVec.y * 3 + sendVec.z * 9 + 13;
       for ( int rank = 0; rank < comm.ensemble().size(); ++rank ) {
         int start = _ens_gatherer->Ens_g_RecvStart( rank, ens_buffer_idx );
         int num = _ens_gatherer->Ens_g_RecvCount( rank, ens_buffer_idx, ens_buffer_idx );
         sendbuffer.CopyFrom( g_ensRecvBuf, num, start, sendbuffer.Number() );
       }
     }

    }

    Vec3<int> recvVec( sendVec );
    recvVec[direction] = 0;
    auto& recvbuffer = GetSendBuffer<P>()[ BufferVecToBufferIdx(recvVec) ];

    auto reqs = MPIHelper::null_requests( 2 );
    if ( dest_rank != NEIGHBOR_NULL ) {
      reqs[0] = comm.cartesian().Isend( dest_rank, 0, sendbuffer.PtcData().data(), sendbuffer.Number() );
    }
    if ( src_rank != NEIGHBOR_NULL ) {
      reqs[1] = comm.cartesian().Irecv( src_rank, 0, recvbuffer.PtcData().data() + recvbuffer.Number(), _sizeParticleBuffer );
    }

    // Wait and then get the recv counts
    auto statuses = MPIHelper::waitall( reqs );

    if( src_rank != NEIGHBOR_NULL ) {
      int count_recv = MPIHelper::get_count<decltype(recvbuffer.PtcData()[0])>( statuses[1] );
      // change the number of recvArray
      recvbuffer.SetNum( recvbuffer.Number() + count_recv );
    }

    //zero out sendbuffer and set number to 0
    if ( dest_rank != NEIGHBOR_NULL ) {
      sendbuffer.Erase( 0, sendbuffer.Number() );
      sendbuffer.SetNum(0);
    }
}

// copy newly arrived particles in the CENTER buffer to particle array
template < typename P >
void ParticleCommunicator::ImportRecvedParticles(Particles<P> &particles, const MPIEnsembleCommunicator& ensemble){
  if ( ensemble.is_root() ) {
    auto& center_buffer = GetSendBuffer<P>()[static_cast<int>(CENTER)];
    auto num_center = center_buffer.Number();

    if ( ensemble.size() > 1 && _ens_rank_import != ensemble.rank() ) {
      ensemble.send( _ens_rank_import, 147, center_buffer.PtcData().data(), num_center );
    } else {
      particles.CopyFrom( center_buffer, num_center, 0, particles.Number() );
    }

    center_buffer.Erase( 0, num_center );
    center_buffer.SetNum(0);

  } else if ( _ens_rank_import == ensemble.rank() ) {
    MPI_Status status = ensemble.probe( ensemble.root(), 147 );
    int num_recv = MPIHelper::get_count<decltype(particles.PtcData()[0])>( status );
    ensemble.recv( ensemble.root(), 147, particles.PtcData().data() + particles.Number(), num_recv );
    particles.SetNum( particles.Number() + num_recv );
  }

}

template < typename P >
void AdjustCells(Particles<P> &particles, const Grid& grid){
  ;
  for( unsigned int i=0; i<particles.Number(); i++ ){
      if( particles.IsEmpty(i) ) continue;

      //assume all nodes have the same grid.dims
      int cell = particles.PtcData()[i].cell;
      int C1 = ( grid.getC1(cell)-grid.guard[0]+grid.reducedDim(0) )%(grid.reducedDim(0)) + grid.guard[0];
      int C2 = grid.dimension > 1 ?
            ( grid.getC2(cell)-grid.guard[1]+grid.reducedDim(1) )%(grid.reducedDim(1)) + grid.guard[1] : 0;
      int C3 = grid.dimension > 2 ?
            ( grid.getC3(cell)-grid.guard[2]+grid.reducedDim(2) )%(grid.reducedDim(2)) + grid.guard[2] : 0;

      particles.PtcData()[i].cell = grid.getIdx(C1,C2,C3);
  }
}

int CountSentParticles( const ParticleSorter& sorter, const std::array<int,3>& neighbor_left, const std::array<int,3>& neighbor_right ) {
  int totalSend = 0;
  const auto& numInTile = sorter.getNumInTile();

  int min0 = (neighbor_left[0] == NEIGHBOR_NULL)? 0 : -1;
  int max0 = (neighbor_right[0] == NEIGHBOR_NULL)? 0 : 1;

  int min1 = (neighbor_left[1] == NEIGHBOR_NULL)? 0 : -1;
  int max1 = (neighbor_right[1] == NEIGHBOR_NULL)? 0 : 1;

  int min2 = (neighbor_left[2] == NEIGHBOR_NULL)? 0 : -1;
  int max2 = (neighbor_right[2] == NEIGHBOR_NULL)? 0 : 1;

  for ( int k = min2; k <= max2; ++k ) {
    for ( int j = min1; j <= max1; ++j ) {
      for ( int i = min0; i <= max0; ++i ) {
        if ( 0 == i && 0 == j && 0 == k ) continue;
        int tile = BufferVecToTile( { i, j, k } );
        totalSend += numInTile[tile];
      }
    }
  }

  return totalSend;
}

template < typename P >
void ParticleCommunicator::SendRecvParticles(int timestep, Particles<P> &particles, const ParticleSorter& sorter, const MPICommunicator& comm, const AperParams& params) {
  const auto& cartesian = comm.cartesian();
  const auto& ensemble = comm.ensemble();
  auto& contents = InfoCollector::Instance().contents;
  const auto& grid = params.grid;
  int ptcType_idx = static_cast<int>(particles.Attributes().ptcType);
  const auto& neighbor_left = params.ens_specs.neighbor_left;
  const auto& neighbor_right = params.ens_specs.neighbor_right;

  // register number of particles to be sent away before actual communication takes place
  int totalsend = CountSentParticles( sorter, neighbor_right, neighbor_right );

  if ( ensemble.size() > 1 ) {
    ensemble.barrier(); // for better timing on ens comm
    auto t0 = high_resolution_clock::now();
    _ens_gatherer->EnsGatherParticles( particles, ensemble, sorter );
    if ( comm.is_primary() ) {
      // record number of gathered particles
      contents.num_ens_gatherPtc[ptcType_idx] = _ens_gatherer->CountTotalGatheredPtcs();
    }
    auto t1 = high_resolution_clock::now();
    auto dur = duration_cast<clock_cycle>( t1 - t0 );
    contents.t_ens_gatherPtc[ptcType_idx] = dur.count();

    _ens_rank_import = timestep % ensemble.size();
  }

  if ( comm.is_primary() ) {
    //front and back in 3D
    if( grid.dimension > 2 ){
      for( int k=-1; k<=1; k+=2 ){
        for( int j=-1; j<=1; j++ ){
          for( int i=-1; i<=1; i++ ){
            Vec3<int> sendVec(i,j,k);
            SendRecvDirectional(particles, sorter, comm, sendVec, 2, neighbor_left[2], neighbor_right[2]);
            cartesian.barrier();
          }
        }
      }
    }

    //up and down
    if ( grid.dimension > 1 ) {
      for( int j=-1; j<=1; j+=2 ){
        for( int i=-1; i<=1; i++ ){
          Vec3<int> sendVec(i,j,0);
          SendRecvDirectional(particles, sorter, comm, sendVec, 1, neighbor_left[1], neighbor_right[1]);
          cartesian.barrier();
        }
      }
    }

    //left and right
    for( int i=-1; i<=1; i+=2 ){
      Vec3<int> sendVec(i,0,0);
      SendRecvDirectional(particles, sorter, comm, sendVec, 0, neighbor_left[0], neighbor_right[0]);
      cartesian.barrier();
    }

    //Recalibrate and copy particles from recv buffer( i.e. _bufferSendParticles[CENTER])
    //to the particle array
    AdjustCells( GetSendBuffer<P>()[static_cast<int>(CENTER)], grid );
  }

  //remove sent-away particles
  const std::vector<int>& numInTile = sorter.getNumInTile();
  const std::vector<int>& scanTile = sorter.getScanTile();
  int tile_center = BufferVecToTile( { 0, 0, 0 } );
  // the -1 is because of no need to erase empty particles
  for( unsigned int i = 0; i < sorter.getNumTiles() - 1; i++ ) {
    if( tile_center == i ) continue; // 13 = bulk region
    particles.Erase( scanTile[i], numInTile[i] );
  }

  // total number recved is just the Number of the SendBuffer CENTER
  int totalrecv= GetSendBuffer<P>()[static_cast<int>(CENTER)].Number();

  auto t_import0 = high_resolution_clock::now();
  ImportRecvedParticles( particles, comm.ensemble() );
  auto t_import1 = high_resolution_clock::now();
  auto dur = duration_cast<clock_cycle>( t_import1 - t_import0 );
  contents.t_importPtcs[ptcType_idx] = dur.count();

  Logger::print( 0, "------ RankCart ", cartesian.rank(), ": Total # of sent", particles.NameStr(), "= ", totalsend, ", Total # of received", particles.NameStr(), "= ", totalrecv, " ------" );
  contents.num_send[ptcType_idx] = totalsend;
  contents.num_recv[ptcType_idx] = totalrecv;
}

void ParticleCommunicator::Reset( int ens_size ) {
  _ens_gatherer->AdjustBufferSizes( ens_size );
}

// cannot partial specialize a function template, so temporarily remove all blocking sendrecv
// NOTE only need to instantiate public members
#define INSTANTIATE_PTC_COMM(P)    \
  template void ParticleCommunicator::SendRecvParticles<P>(int timestep, Particles<P>& particles, const ParticleSorter& sorter, const MPICommunicator& comm, const AperParams& params ); \

INSTANTIATE_PTC_COMM(QCarrier);
INSTANTIATE_PTC_COMM(QNeutral);
