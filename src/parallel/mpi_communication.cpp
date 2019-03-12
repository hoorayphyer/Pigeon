#include "parallel/mpi_communication.hpp"
#include <type_traits>
#include <experimental/type_traits> // for is_detected
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
  template <typename T>
  using is_container_t = decltype( std::declval<T>().data()[  std::declval<T>().size() ] );

  template <typename T>
  constexpr bool is_container() {
    return std::experimental::is_detected< is_container_t, T >::value;
  }

  template <typename T>
  inline auto decay_buf( T&& buf ) {
    if constexpr ( is_container<T>() )
                   return std::make_tuple( buf.data(), buf.size(), datatype(buf.data()) );
    else
      return std::make_tuple( &buf, 1, datatype(&buf) );
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

  constexpr auto mpi_op( ReduceOp op ) {
    if ( ReduceOp::SUM == op ) return MPI_SUM;
    if ( ReduceOp::MAX == op ) return MPI_MAX;
    if ( ReduceOp::MAXLOC == op ) return MPI_MAXLOC;
  }

  template < bool In_Place, typename T >
  std::tuple< const void*, void* > find_out_reduce_buffers( void* buf, const T& buffer, reduce_return_t<T>& result, bool is_root ) {
    const void* send_buf = nullptr;
    void* recv_buf = nullptr;

    if ( !is_root ) {
      send_buf = buf;
    } else {
      if constexpr( In_Place ) {
          send_buf = MPI_IN_PLACE;
          recv_buf = buf;
        } else {
        send_buf = buf;
        result.emplace( buffer ); // copy construct the recv_buf
        recv_buf = std::get<0>( decay_buf( *result ) );
      }
    }

    return std::make_tuple( send_buf, recv_buf );
  }

  template < typename Comm >
  template < bool In_Place, ReduceOp op, typename T>
  reduce_return_t<T> Collective_Comm<Comm>::reduce( T& buffer, int root ) const {
    reduce_return_t<T> result;

    auto[ buf, count, datatype ] = decay_buf( std::forward<T>(buffer) );
    auto[ send_buf, recv_buf ] = find_out_reduce_buffers<In_Place>( buf, buffer, result, _comm().rank() == root );

    MPI_Reduce( send_buf, recv_buf, count, datatype, mpi_op(op), root, _comm() );

    return result;
  }

  template < typename Comm >
  template < bool In_Place, ReduceOp op, typename T>
  std::tuple<Request, reduce_return_t<T> > Collective_Comm<Comm>::Ireduce( T& buffer, int root ) const {
    reduce_return_t<T> result;

    auto[ buf, count, datatype ] = decay_buf( std::forward<T>(buffer) );
    auto[ send_buf, recv_buf ] = find_out_reduce_buffers<In_Place>( buf, buffer, result, _comm().rank() == root );

    Request req;
    MPI_Ireduce( send_buf, recv_buf, count, datatype, mpi_op(op), root, _comm(), req );

    return std::make_tuple( req, result );
  }

  template < typename Comm>
  template < typename T >
  void Collective_Comm<Comm>::broadcast( T& buffer, int root ) const {
    auto[ buf, count, datatype ] = decay_buf(buffer);
    MPI_Bcast( buf, count, datatype, root, _comm() );
  }


  template < typename Comm>
  template < typename T >
  Request Collective_Comm<Comm>::Ibroadcast( T& buffer, int root ) const {
    Request req;
    auto[ buf, count, datatype ] = decay_buf(buffer);
    MPI_Ibcast( buf, count, datatype, root, _comm(), req );
    return req;
  }
}
