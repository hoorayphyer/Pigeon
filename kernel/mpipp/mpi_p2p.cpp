#include "mpipp/mpi++.hpp"
// TODOL Edge case: send recv to self. Both blocking and nonblocking are for communication between different processes

namespace mpi {
  template < typename Comm >
  template < typename T >
  int P2P_Comm<Comm>::probe( int source_rank, int tag, const T* ) const {
    MPI_Status s;
    MPI_Probe( source_rank, tag, _comm(), &s );
    int count = 0;
    MPI_Get_count( &s, datatype((const T*)0), &count );
    return count;
  }

  template < bool nonblocking >
  constexpr auto pick_send ( SendMode mode ) noexcept {
    switch ( mode ) {
    case SendMode::BUF :
      if constexpr ( nonblocking ) return MPI_Ibsend;
      else return MPI_Bsend;
    case SendMode::SYN :
      if constexpr ( nonblocking ) return MPI_Issend;
      else return MPI_Ssend;
    case SendMode::RDY :
      if constexpr ( nonblocking ) return MPI_Irsend;
      else return MPI_Rsend;
    default :
      if constexpr ( nonblocking ) return MPI_Isend;
      else return MPI_Send;
    }
  }

  // TODOL consider add error handling for send
  template < typename Comm >
  template < typename T >
  void P2P_Comm<Comm>::send(int dest_rank, int tag, const T* send_buf, int send_count, SendMode mode ) const {
    (*(pick_send<false>(mode))) ( send_buf, send_count, datatype<T>(), dest_rank, tag, _comm() );
  }


  template < typename Comm >
  template < typename T >
  Request P2P_Comm<Comm>::Isend( int dest_rank, int tag, const T* send_buf, int send_count, SendMode mode ) const {
    Request req;
    req.reset( new MPI_Request );
    (*(pick_send<true>(mode))) ( send_buf, send_count, datatype<T>(), dest_rank, tag, _comm(), req );

    return req;
  }

  template < typename Comm >
  template < typename T >
  int P2P_Comm<Comm>::recv( int source_rank, int tag, T* recv_buf, int recv_count_max ) const {
    MPI_Status status;
    MPI_Recv( recv_buf, recv_count_max, datatype<T>(), source_rank, tag, _comm(), &status);

    int recv_count = 0;
    MPI_Get_count( &status, datatype<T>(), &recv_count );

    return recv_count;
  }

  template < typename Comm >
  template < typename T >
  Request P2P_Comm<Comm>::Irecv(int source_rank, int tag, T* recv_buf, int recv_count_max ) const {
    Request req;
    req.reset( new MPI_Request );
    MPI_Irecv( recv_buf, recv_count_max, datatype<T>(), source_rank, tag, _comm(), req );
    return req;
  }

}

// TODOL together with mpi_datatype, use typelist to improve
// instantiation
namespace mpi {
#define INSTANTIATE_MPI_P2P_FOR_COMM(_COMM_, _TYPE_)                    \
  template int P2P_Comm<_COMM_>::probe( int, int, const _TYPE_ * ) const; \
  template void P2P_Comm<_COMM_>::send( int, int, const _TYPE_ *, int, SendMode ) const; \
  template Request P2P_Comm<_COMM_>::Isend( int, int, const _TYPE_ *, int, SendMode ) const; \
  template int P2P_Comm<_COMM_>::recv( int, int, _TYPE_ *, int ) const; \
  template Request P2P_Comm<_COMM_>::Irecv( int, int, _TYPE_ *, int ) const

#define INSTANTIATE_MPI_P2P(_TYPE_)                  \
  INSTANTIATE_MPI_P2P_FOR_COMM(Comm, _TYPE_);         \
  INSTANTIATE_MPI_P2P_FOR_COMM(InterComm, _TYPE_)

  INSTANTIATE_MPI_P2P(char);
  INSTANTIATE_MPI_P2P(int);
  INSTANTIATE_MPI_P2P(unsigned int);
  INSTANTIATE_MPI_P2P(float);
  INSTANTIATE_MPI_P2P(double);
  INSTANTIATE_MPI_P2P(long double);
}

#include "mpi_cparticle.cpp"
namespace mpi {
  INSTANTIATE_MPI_P2P(cPtc_t);
}
