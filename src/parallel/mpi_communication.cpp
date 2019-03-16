#include "parallel/mpi_communication.hpp"
#include <type_traits>

namespace mpi {
  void request_free ( MPI_Request* p ) {
    if ( p && *p != MPI_REQUEST_NULL )
      MPI_Request_free(p);
  }

  MPI_Request request_null() {
    return MPI_REQUEST_NULL;
  }

  void wait( Request& request ) {
    // MPI_Status status;
    // MPI_Wait( &request, &status );
    // return status;
    MPI_Wait( request, MPI_STATUS_IGNORE );
  }

  void waitall( std::vector<Request>& requests ) {
    MPI_Request* p = new MPI_Request [ requests.size() ];
    for ( int i = 0; i < requests.size(); ++i )
      p[i] = requests[i];
    // std::vector<MPI_Status> statuses(requests.size());
    // MPI_Waitall( requests.size(), requests.data(), statuses.data() );
    // return statuses;
    MPI_Waitall( requests.size(), p, MPI_STATUSES_IGNORE );
    delete [] p;
  }

  void cancel( Request& req ) {
    if ( static_cast<MPI_Request>(req) != MPI_REQUEST_NULL ) // this is necessary
      MPI_Cancel(req);
  }

  void cancelall( std::vector<Request>& reqs ) {
    for ( auto& req : reqs )
      cancel(req);
  }
}

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
    MPI_Irecv( recv_buf, recv_count_max, datatype<T>(), source_rank, tag, _comm(), req );
    return req;
  }

}

namespace mpi {
  template < typename Comm ,bool Inter >
  void Collective_Comm<Comm, Inter>::barrier() const {
    MPI_Barrier( _comm() );
  }

  constexpr auto mpi_op( by op ) {
    if ( by::SUM == op ) return MPI_SUM;
    if ( by::MAX == op ) return MPI_MAX;
    if ( by::MAXLOC == op ) return MPI_MAXLOC;
  }

  template < typename Comm, bool Inter >
  template < by Op, bool In_Place, typename T >
  std::optional<std::vector<T>>
  Collective_Comm<Comm, Inter>::reduce( T* buffer, int count, int root ) const {
    static_assert( !Inter );
    std::optional<std::vector<T>> result;

    const void* send_buf = nullptr;
    void* recv_buf = nullptr;

    if ( _comm().rank() != root ) {
      send_buf = buffer;
    } else {
      if constexpr( In_Place ) {
          send_buf = MPI_IN_PLACE;
          recv_buf = buffer;
        } else {
        send_buf = buffer;
        result.emplace(); // copy construct the recv_buf
        (*result).resize(count);
        (*result).shrink_to_fit();
        recv_buf = (*result).data();
      }
    }

    MPI_Reduce( send_buf, recv_buf, count, datatype, mpi_op(Op), root, _comm() );

    return result;
  }

  // template < typename Comm, bool Inter >
  // template < by Op, bool In_Place, typename T >
  // std::tuple< Request, std::optional<std::vector<T> > >
  // Collective_Comm<Comm,Inter>::Ireduce( T& buffer, int root ) const {
  //   static_assert( !Inter );
  //   std::optional<std::vector<T>> result;

  //   const void* send_buf = nullptr;
  //   void* recv_buf = nullptr;

  //   if ( _comm().rank() != root ) {
  //     send_buf = buf;
  //   } else {
  //     if constexpr( In_Place ) {
  //         send_buf = MPI_IN_PLACE;
  //         recv_buf = buffer;
  //       } else {
  //       send_buf = buffer;
  //       result.emplace(); // copy construct the recv_buf
  //       (*result).resize(count);
  //       (*result).shrink_to_fit();
  //       recv_buf = (*result).data();
  //     }
  //   }

  //   Request req;
  //   MPI_Ireduce( send_buf, recv_buf, count, datatype, mpi_op(op), root, _comm(), req );

  //   return std::make_tuple( req, result );
  // }

  template < typename Comm, bool Inter>
  template < typename T >
  void Collective_Comm<Comm, Inter>::broadcast( int root, T* buffer, int count ) const {
    static_assert( !Inter );
    MPI_Bcast( buffer, count, datatype<T>(), root, _comm() );
  }


  template < typename Comm, bool Inter >
  template < typename T >
  Request Collective_Comm<Comm, Inter>::Ibroadcast( int root, T* buffer, int count ) const {
    static_assert( !Inter );
    Request req;
    MPI_Ibcast( buffer, count, datatype<T>(), root, _comm(), req );
    return req;
  }


  template < typename Comm, bool Inter >
  template < typename T >
  std::vector<T> Collective_Comm<Comm, Inter>::allgather( const T* send_buf, int send_count ) const {
    std::vector<T> recv( send_count * _comm().size() ); // CLAIM: also works with intercomm
    recv.shrink_to_fit();
    MPI_Allgather( send_buf, send_count, datatype<T>(), recv.data(), recv.size(), datatype<T>(), _comm() );
    return recv;

  }

  template < typename Comm, bool Inter>
  template < typename T >
  void Collective_Comm<Comm, Inter>::scatter( int root, T* buffer, int count ) const {
    void* sendbuf = nullptr;
    void* recvbuf = nullptr;
    int sendcount = 0, recvcount = 0;
    if constexpr ( Inter ) {
        if ( MPI_ROOT == root ) {
          sendbuf = buffer;
          sendcount = count;
        }
        else if ( MPI_PROC_NULL == root );
        else {
          recvbuf = buffer;
          recvcount = count;
        }
      } else {
      if ( _comm().rank() == root ) {
        sendbuf = buffer;
        sendcount = count;
        recvbuf = MPI_IN_PLACE;
      } else {
        recvbuf = buffer;
        recvcount = count;
      }
    }

    MPI_Scatter(sendbuf, sendcount, datatype<T>(), recvbuf, recvcount, datatype<T>(), root, _comm() );
  }
}
