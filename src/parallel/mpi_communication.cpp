#include "parallel/mpi_communication.hpp"
#include <type_traits>
#include <mpi.h>

namespace mpi {
  template <typename Type>
  MPI_Datatype datatype(Type*) noexcept {
    using T = std::remove_const_t<Type>;
    if constexpr ( std::is_same_v<T, char> ) return MPI_CHAR;
    else if ( std::is_same_v<T, short> ) return MPI_SHORT;
    else if ( std::is_same_v<T, int> ) return MPI_INT;
    else if ( std::is_same_v<T, long> ) return MPI_LONG;

    else if ( std::is_same_v<T, unsigned char> ) return MPI_UNSIGNED_CHAR;
    else if ( std::is_same_v<T, unsigned short> ) return MPI_UNSIGNED_SHORT;
    else if ( std::is_same_v<T, unsigned int> ) return MPI_UNSIGNED;
    else if ( std::is_same_v<T, unsigned long> ) return MPI_UNSIGNED_LONG;

    else if ( std::is_same_v<T, float> ) return MPI_FLOAT;
    else if ( std::is_same_v<T, double> ) return MPI_DOUBLE;

    else if ( std::is_same_v<T, long double> ) return MPI_LONG_DOUBLE;
    else if ( std::is_same_v<T, bool> ) return MPI_CXX_BOOL;
    else return MPI_DATATYPE_NULL;
  }
}

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
  template < typename Comm >
  void Collective_Comm<Comm>::barrier() const {
    MPI_Barrier( _comm() );
  }

  constexpr auto mpi_op( by op ) {
    if ( by::SUM == op ) return MPI_SUM;
    if ( by::MAX == op ) return MPI_MAX;
    if ( by::MAXLOC == op ) return MPI_MAXLOC;
  }

  template < typename Comm >
  template < by Op, bool In_Place, typename T >
  std::optional<std::vector<T>>
  Collective_Comm<Comm>::reduce( T* buffer, int count, int root ) const {
    std::optional<std::vector<T>> result;

    const void* send_buf = nullptr;
    void* recv_buf = nullptr;

    if ( _comm().rank() != root ) {
      send_buf = buf;
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

    MPI_Reduce( send_buf, recv_buf, count, datatype, mpi_op(op), root, _comm() );

    return result;
  }

  template < typename Comm >
  template < by Op, bool In_Place, typename T >
  std::tuple<Request, reduce_return_t<T> > Collective_Comm<Comm>::Ireduce( T& buffer, int root ) const {
    std::optional<std::vector<T>> result;

    const void* send_buf = nullptr;
    void* recv_buf = nullptr;

    if ( _comm().rank() != root ) {
      send_buf = buf;
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

    Request req;
    MPI_Ireduce( send_buf, recv_buf, count, datatype, mpi_op(op), root, _comm(), req );

    return std::make_tuple( req, result );
  }

  template < bool In_Place, ReduceOp Op, typename T >
  std::tuple< Request, std::optional<std::vector<T> > >
  Ireduce( T* buffer, int count, int root ) const;

  template < typename Comm>
  template < typename T >
  void Collective_Comm<Comm>::broadcast( int root, T* buffer, int count ) const {
    MPI_Bcast( buffer, count, datatype<T>(), root, _comm() );
  }


  template < typename Comm>
  template < typename T >
  Request Collective_Comm<Comm>::Ibroadcast( int root, T* buffer, int count ) const {
    Request req;
    MPI_Ibcast( buffer, count, datatype<T>(), root, _comm(), req );
    return req;
  }


  template < typename Comm>
  template < typename T >
  std::vector<T> Collective_Comm<Comm>::allgather( const T* send_buf, int send_count ) const {
    std::vector<T> recv( send_count * _comm().size() ); // CLAIM: also works with intercomm
    recv.shrink_to_fit();
    MPI_Allgather( send_buf, send_count, datatype<T>(), recv.data(), recv.size(), datatype<T>(), _comm() );
    return recv;

  }
}
