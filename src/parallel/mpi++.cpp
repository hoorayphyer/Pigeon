#include "parallel/mpi++.hpp"
#include <mpi.h>
#include <algorithm>
#include <iterator> // back_inserter
#include <tuple>
#include <numeric> // std::inner_product

// helpers
namespace mpi {
  template< typename MPI_Handle_t>
  struct Raw;

  template<>
  struct Raw<MPI_Request> {
    constexpr static auto null_hdl = MPI_REQUEST_NULL;
    constexpr static void free( MPI_Request* p ) {
      if ( !p && *p != null_hdl ) MPI_Request_free(p);
    }
  };

  // template<>
  // struct Raw<MPI_Group> {
  //   constexpr static auto null_hdl = MPI_GROUP_NULL;
  //   constexpr static void free(MPI_Group* p) {
  //     if ( !p && *p != null_hdl ) MPI_Group_free(p);
  //   }
  // };

  template<>
  struct Raw<MPI_Comm> {
    constexpr static auto null_hdl = MPI_COMM_NULL;
    constexpr static void free(MPI_Comm* p) {
      // attempt to free MPI_COMM_WORLD is forbidden
      if ( !p && *p != null_hdl && *p != MPI_COMM_WORLD ) MPI_Comm_free(p);
    }
  };

}

// global definitions
namespace mpi {
  Comm world;

  void initialize() {
    int is_initialized = 0;
    MPI_Initialized(&is_initialized);

    if (!is_initialized)
      MPI_Init(NULL, NULL);

    MPI_Comm comm_world = MPI_COMM_WORLD;
    world.hdl = Handle(&comm_world, Raw<MPI_Comm>::free );
  }

  void finalize() {
    int is_finalized = 0;
    MPI_Finalized(&is_finalized);

    if (!is_finalized) MPI_Finalize();
  }
}

// mpi::Group
namespace mpi {
  Group operator^ ( Group a, Group b ) {
    Group result;
    std::set_intersection( a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(result) );
    return result;
  }

  Group operator+ ( Group a, Group b ) {
    Group result;
    std::set_union( a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(result) );
    return result;
  }

  Group operator- ( Group a, Group b ) {
    Group result;
    std::set_difference( a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(result) );
    return result;
  }
}

// mpi::Comm
namespace mpi {
  Comm::Comm( Group group ) {
    MPI_Group grp_world, grp_this;
    MPI_Comm_group(MPI_COMM_WORLD, &grp_world);
    MPI_Group_incl( grp_world, group.size(), group.data(), &grp_this );

    MPI_Comm comm;
    MPI_Comm_create_group(MPI_COMM_WORLD, grp_this, 0, &comm );
    hdl = Handle( &comm, Raw<MPI_Comm>::free );

    MPI_Group_free(&grp_world);
    MPI_Group_free(&grp_this);
  }

  int Comm::rank() const {
    int rank;
    MPI_Comm_rank( hdl, &rank );
    return rank;
  }

  int Comm::size() const {
    int size;
    MPI_Comm_size( hdl, &size );
    return size;
  }

  std::vector<int> Comm::to_world( std::vector<int> ranks ) const {
    auto world_ranks {ranks};
    MPI_Group grp_world, grp_this;
    MPI_Comm_group(MPI_COMM_WORLD, &grp_world);
    MPI_Comm_group(hdl, &grp_this);
    MPI_Group_translate_ranks( grp_this, ranks.size(), ranks.data(), grp_world, world_ranks.data() );
    MPI_Group_free( &grp_world );
    MPI_Group_free( &grp_this );
    return world_ranks;
  }

}

// mpi::topo
namespace mpi {
  namespace topo {
    void cartesianize( Comm& comm, std::vector<int> dims, std::vector<bool> periodic ) {
      const int ndims = dims.size() < periodic.size() ? dims.size() : periodic.size();

      std::vector<int> periods(ndims);
      for ( int i = 0; i < ndims; ++i )
        periods[i] = static_cast<int>( periodic[i] );

      MPI_Comm comm_cart;
      MPI_Cart_create( comm.hdl, ndims, dims.data(), periods.data(), true, &comm_cart );
      comm.hdl.reset( &comm_cart );
    }

    namespace cart {
      std::tuple<std::vector<int>, std::vector<bool>, std::vector<int>> inquire ( const Comm& comm ) {
        int ndims;
        MPI_Cartdim_get( comm.hdl, &ndims );
        std::vector<int> dims(ndims);
        std::vector<int> periods(ndims);
        std::vector<int> coords(ndims);
        MPI_Cart_get( comm.hdl, ndims, dims.data(), periods.data(), coords.data() );

        std::vector<bool> periodic( ndims );
        for ( int i = 0; i < ndims; ++i )
          periodic[i] = static_cast<bool>( periods[i] );

        return std::make_tuple( dims, periodic, coords );
      }

      std::vector<int> dims ( const Comm& comm ) {
        return std::get<0>( inquire( comm ) );
      }

      std::vector<bool> periodic ( const Comm& comm ) {
        auto periods = std::get<1>( inquire( comm ) );
        std::vector<bool> result(periods.size()); // convert vector<int> to vector<bool>
        std::copy( periods.begin(), periods.end(), result.begin() );
        return result;
      }

      std::vector<int> coords ( const Comm& comm ) {
        return std::get<2>( inquire( comm ) );
      }

      int linear_coord( const Comm& comm ) {
        auto[strides, to_be_ignored, coords] = inquire(comm);
        // do exclusive multiplication on stride to get real strides
        strides.back() = std::accumulate( strides.begin(), strides.end() - 1, 1, [](auto a, auto b){return a*b;} );
        for ( int i = strides.size() - 2; i >= 0; --i )
          strides[i] = strides[i+1] / strides[i];

        return std::inner_product( coords.begin(), coords.end(), strides.begin(), 0 );
      }


    }
  }
}


#include <experimental/type_traits> // for is_detected
// mpi::Comm communication functions
namespace mpi {
  template <typename T>
  using is_container_t = decltype( std::declval<T>().data()[  std::declval<T>().size() ] );

  template <typename T>
  constexpr bool is_container() {
    return std::experimental::is_detected< is_container_t, T >::value;
  }

  template <typename T>
  constexpr MPI_Datatype mpi_datatype() {
    // TODO
    // using typename std::remove_const<typename std::remove_reference<T>::type>::type
    return MPI_INT;
  }



  template <typename T>
  inline auto decode_buf( T&& buf ) {
    if constexpr ( is_container<T>(buf) )
      return std::make_tuple( buf.data(), buf.size(), mpi_datatype<decltype(buf[0])>() );
    else
      return std::make_tuple( &buf, 1, mpi_datatype<decltype(buf)>() );
  }


  void Comm::barrier() const {
    MPI_Barrier( hdl );
  }


  template < bool nonblocking >
  constexpr auto pick_send ( SendMode mode ) {
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
  template <typename T>
  void Comm::send(int dest_rank, const T& send_buf, int tag, SendMode mode ) const {
    auto[ buf, count, datatype ] = decode_buf(send_buf);
    *(pick_send<false>(mode)) ( buf, count, datatype, dest_rank, tag, hdl );
  }


  template <typename T>
  Request Comm::Isend( int dest_rank, const T& send_buf, int tag, SendMode mode ) const {
    Request req;
    req.hdl = Handle( new MPI_Request, Raw<MPI_Request>::free );

    auto[ buf, count, datatype ] = decode_buf(send_buf);
    *(pick_send<true>(mode)) ( buf, count, datatype, dest_rank, tag, hdl, &req.hdl );

    return req;
  }

  int probe_size( int source_rank, int tag, const MPI_Comm& comm, MPI_Datatype datatype  ) {
    MPI_Status s;
    MPI_Probe( source_rank, tag, comm, &s );
    int count = 0;
    MPI_Get_count( &s, datatype, &count );
    return count;
  }

  template <typename T>
  int Comm::recv( int source_rank, T& recv_buf, int tag, bool resize_buf_with_probe ) const {
    int recv_count = 0;

    auto[ buf, count, datatype ] = decode_buf( recv_buf );

    if constexpr ( is_container<T>() ) {
      if ( resize_buf_with_probe ) {
        count = probe_size( source_rank, tag, hdl, datatype );
        recv_buf.resize(count);
      }
    }

    MPI_Status status;
    MPI_Recv( buf, count, datatype, source_rank, tag, hdl, &status);
    MPI_Get_count( &status, datatype, &recv_count );

    return recv_count;
  }

  template <typename T>
  Request Comm::Irecv(int source_rank, T& recv_buf, int tag, bool resize_buf_with_probe ) const {
    Request req;
    req.hdl = Handle( new MPI_Request, Raw<MPI_Request>::free );

    auto[ buf, count, datatype ] = decode_buf( recv_buf );

    if constexpr ( is_container<T>() ) {
      if ( resize_buf_with_probe ) {
        count = probe_size( source_rank, tag, hdl, datatype );
        recv_buf.resize(count);
      }
    }

    MPI_Irecv( recv_buf.data(), recv_buf.size(), mpi_datatype<decltype(recv_buf[0])>(), source_rank, tag, hdl, &req.hdl );
    return req;
  }

  constexpr auto mpi_op( ReduceOp op ) {
    if ( ReduceOp::SUM == op ) return MPI_SUM;
    if ( ReduceOp::MAX == op ) return MPI_MAX;
    if ( ReduceOp::MAXLOC == op ) return MPI_MAXLOC;
  }

  template < bool In_Place, typename T >
  std::tuple< const void*, void* > find_out_reduce_buffers( void* buf, const T& buffer, Comm::reduce_return_t<T>& result, bool is_root ) {
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
        recv_buf = std::get<0>( decode_buf( *result ) );
      }
    }

    return std::make_tuple( send_buf, recv_buf );
  }

  template < bool In_Place, ReduceOp op, typename T>
  Comm::reduce_return_t<T> Comm::reduce( T& buffer, int root ) const {
    reduce_return_t<T> result;

    auto[ buf, count, datatype ] = decode_buf( std::forward<T>(buffer) );
    auto[ send_buf, recv_buf ] = find_out_reduce_buffers<In_Place>( buf, buffer, result, rank() == root );

    MPI_Reduce( send_buf, recv_buf, count, datatype, mpi_op(op), root, hdl );

    return result;
  }

  template < bool In_Place, ReduceOp op, typename T>
  std::tuple<Request, Comm::reduce_return_t<T> > Comm::Ireduce( T& buffer, int root ) const {
    reduce_return_t<T> result;

    auto[ buf, count, datatype ] = decode_buf( std::forward<T>(buffer) );
    auto[ send_buf, recv_buf ] = find_out_reduce_buffers<In_Place>( buf, buffer, result, rank() == root );

    Request req;
    req.hdl = Handle ( new MPI_Request, Raw<MPI_Request>::free );

    MPI_Ireduce( send_buf, recv_buf, count, datatype, mpi_op(op), root, hdl, &req.hdl );

    return std::make_tuple( req, result );
  }

  template < typename T >
  void Comm::broadcast( T& buffer, int root ) const {
    auto[ buf, count, datatype ] = decode_buf(buffer);
    MPI_Bcast( buf, count, datatype, root, hdl );
  }


  template < typename T >
  Request Comm::Ibroadcast( T& buffer, int root ) const {
    Request req;
    req.hdl = Handle( new MPI_Request, Raw<MPI_Request>::free );

    auto[ buf, count, datatype ] = decode_buf(buffer);
    MPI_Ibcast( buf, count, datatype, root, hdl, &req.hdl );

    return req;
  }
}
