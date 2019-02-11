#include "parallel/mpi_communication.hpp"
#include "./mpi_datatype.cpp"
#include <experimental/type_traits> // for is_detected
#include <mpi.h>

namespace mpi {
  template <typename T>
  using is_container_t = decltype( std::declval<T>().data()[  std::declval<T>().size() ] );

  template <typename T>
  constexpr bool is_container() {
    return std::experimental::is_detected< is_container_t, T >::value;
  }

  template <typename T>
  inline auto decay_buf( T&& buf ) {
    if constexpr ( is_container<T>(buf) )
      return std::make_tuple( buf.data(), buf.size(), mpi_datatype<decltype(buf[0])>() );
    else
      return std::make_tuple( &buf, 1, mpi_datatype<decltype(buf)>() );
  }

}

namespace mpi {
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
  void P2P_Comm<Comm>::send(int dest_rank, const T& send_buf, int tag, SendMode mode ) const {
    auto[ buf, count, datatype ] = decay_buf(send_buf);
    *(pick_send<false>(mode)) ( buf, count, datatype, dest_rank, tag, _comm().hdl );
  }


  template < typename Comm >
  template < typename T >
  Request P2P_Comm<Comm>::Isend( int dest_rank, const T& send_buf, int tag, SendMode mode ) const {
    Request req;
    req.hdl = Handle( new MPI_Request, Raw<MPI_Request>::free );

    auto[ buf, count, datatype ] = decay_buf(send_buf);
    *(pick_send<true>(mode)) ( buf, count, datatype, dest_rank, tag, _comm().hdl, &req.hdl );

    return req;
  }

  int probe_size( int source_rank, int tag, const MPI_Comm& comm, MPI_Datatype datatype ) {
    MPI_Status s;
    MPI_Probe( source_rank, tag, comm, &s );
    int count = 0;
    MPI_Get_count( &s, datatype, &count );
    return count;
  }

  template < typename Comm >
  template < typename T >
  int P2P_Comm<Comm>::recv( int source_rank, T& recv_buf, int tag, bool resize_buf_with_probe ) const {
    int recv_count = 0;

    auto[ buf, count, datatype ] = decay_buf( recv_buf );

    if constexpr ( is_container<T>() ) {
      if ( resize_buf_with_probe ) {
        count = probe_size( source_rank, tag, _comm().hdl, datatype );
        recv_buf.resize(count);
      }
    }

    MPI_Status status;
    MPI_Recv( buf, count, datatype, source_rank, tag, _comm().hdl, &status);
    MPI_Get_count( &status, datatype, &recv_count );

    return recv_count;
  }

  template < typename Comm >
  template < typename T >
  Request PCP_Comm<Comm>::Irecv(int source_rank, T& recv_buf, int tag, bool resize_buf_with_probe ) const {
    Request req;
    req.hdl = Handle( new MPI_Request, Raw<MPI_Request>::free );

    auto[ buf, count, datatype ] = decay_buf( recv_buf );

    if constexpr ( is_container<T>() ) {
      if ( resize_buf_with_probe ) {
        count = probe_size( source_rank, tag, _comm().hdl, datatype );
        recv_buf.resize(count);
      }
    }

    MPI_Irecv( recv_buf.data(), recv_buf.size(), mpi_datatype<decltype(recv_buf[0])>(), source_rank, tag, _comm().hdl, &req.hdl );
    return req;
  }

}

namespace mpi {
  template < typename Comm >
  void Collective_Comm<Comm>::barrier() const {
    MPI_Barrier( _comm().hdl );
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
    auto[ send_buf, recv_buf ] = find_out_reduce_buffers<In_Place>( buf, buffer, result, rank() == root );

    MPI_Reduce( send_buf, recv_buf, count, datatype, mpi_op(op), root, _comm().hdl );

    return result;
  }

  template < typename Comm >
  template < bool In_Place, ReduceOp op, typename T>
  std::tuple<Request, reduce_return_t<T> > Collective_Comm<Comm>::Ireduce( T& buffer, int root ) const {
    reduce_return_t<T> result;

    auto[ buf, count, datatype ] = decay_buf( std::forward<T>(buffer) );
    auto[ send_buf, recv_buf ] = find_out_reduce_buffers<In_Place>( buf, buffer, result, rank() == root );

    Request req;
    req.hdl = Handle ( new MPI_Request, Raw<MPI_Request>::free );

    MPI_Ireduce( send_buf, recv_buf, count, datatype, mpi_op(op), root, _comm().hdl, &req.hdl );

    return std::make_tuple( req, result );
  }

  template < typename Comm>
  template < typename T >
  void Collective_Comm<Comm>::broadcast( T& buffer, int root ) const {
    auto[ buf, count, datatype ] = decay_buf(buffer);
    MPI_Bcast( buf, count, datatype, root, _comm().hdl );
  }


  template < typename Comm>
  template < typename T >
  Request Collective_Comm<Comm>::Ibroadcast( T& buffer, int root ) const {
    Request req;
    req.hdl = Handle( new MPI_Request, Raw<MPI_Request>::free );

    auto[ buf, count, datatype ] = decay_buf(buffer);
    MPI_Ibcast( buf, count, datatype, root, _comm().hdl, &req.hdl );

    return req;
  }
}
