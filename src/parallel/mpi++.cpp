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

  inline std::shared_ptr<void>& impl_cast( Handle& hdl ) { return hdl._ptr; }
  inline const std::shared_ptr<void>& impl_cast( const Handle& hdl ) { return hdl._ptr; }

  template <typename T, typename Handle_t = Handle>
  inline decltype(auto) as( Handle_t&& hdl ) {
    return *std::static_pointer_cast<T>( impl_cast( std::forward<Handle_t>(hdl) ) );
  }

}

// global definitions
namespace mpi {
  Comm world;
  Handle NULL_HANDLE;

  bool Handle::operator==( Handle other ) const  {
    return _ptr.get() == other._ptr.get(); // TODO check this
  }


  void initialize() {
    int is_initialized = 0;
    MPI_Initialized(&is_initialized);

    if (!is_initialized)
      MPI_Init(NULL, NULL);

    MPI_Comm comm_world = MPI_COMM_WORLD;
    impl_cast(world.hdl).reset( &comm_world, Raw<MPI_Comm>::free );
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
    impl_cast(hdl).reset( &comm, Raw<MPI_Comm>::free );

    MPI_Group_free(&grp_world);
    MPI_Group_free(&grp_this);
  }

  int Comm::rank() const {
    int rank;
    MPI_Comm_rank( as<MPI_Comm>(hdl), &rank );
    return rank;
  }

  int Comm::size() const {
    int size;
    MPI_Comm_size( as<MPI_Comm>(hdl), &size );
    return size;
  }

  std::vector<int> Comm::to_world( std::vector<int> ranks ) const {
    auto world_ranks {ranks};
    MPI_Group grp_world, grp_this;
    MPI_Comm_group(MPI_COMM_WORLD, &grp_world);
    MPI_Comm_group(as<MPI_Comm>(hdl), &grp_this);
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
      MPI_Cart_create( as<MPI_Comm>(comm.hdl), ndims, dims.data(), periods.data(), true, &comm_cart );
      impl_cast(comm.hdl).reset(&comm_cart);
    }

    namespace cart {
      std::tuple<std::vector<int>, std::vector<bool>, std::vector<int>> inquire ( const Comm& comm ) {
        int ndims;
        MPI_Cartdim_get( as<MPI_Comm>(comm.hdl), &ndims );
        std::vector<int> dims(ndims);
        std::vector<int> periods(ndims);
        std::vector<int> coords(ndims);
        MPI_Cart_get( as<MPI_Comm>(comm.hdl), ndims, dims.data(), periods.data(), coords.data() );

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
    // using typename std::remove_const<typename std::remove_reference<T>::type>::type
    return MPI_INT;
  }


  void Comm::barrier() const {
    MPI_Barrier( as<MPI_Comm>(hdl) );
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

  // TODO consider add error handling for send
  template <typename T>
  void Comm::send(int dest_rank, const T& send_buf, int tag, SendMode mode ) const {
    if constexpr ( is_container<T>() )
      *(pick_send<false>(mode)) (send_buf.data(), send_buf.size(), mpi_datatype<decltype(send_buf[0])>(), dest_rank, tag, as<MPI_Comm>(hdl));
    else
      *(pick_send<false>(mode)) (&send_buf, 1, mpi_datatype<decltype(send_buf)>(), dest_rank, tag, as<MPI_Comm>(hdl));
  }


  template <typename T>
  Request Comm::Isend( int dest_rank, const T& send_buf, int tag, SendMode mode ) const {
    Request req;
    MPI_Request mpi_req = MPI_REQUEST_NULL;

    if constexpr ( is_container<T>() )
      *(pick_send<true>(mode))(send_buf.data(), send_buf.size(), mpi_datatype<decltype(send_buf[0])>(), dest_rank, tag, as<MPI_Comm>(hdl), &mpi_req);
    else
      *(pick_send<true>(mode))(&send_buf, 1, mpi_datatype<decltype(send_buf[0])>(), dest_rank, tag, as<MPI_Comm>(hdl), &mpi_req);

    impl_cast(req.hdl).reset( &mpi_req, Raw<MPI_Request>::free );

    return req;
  }

  template <typename T>
  int Comm::recv( int source_rank, T& recv_buf, int tag, bool resize_buf_with_probe ) const {
    int recv_count = 0;

    MPI_Status status;
    if constexpr ( is_container<T>() ) {
      if ( resize_buf_with_probe ) {
        MPI_Status s;
        MPI_Probe( source_rank, tag, as<MPI_Comm>(hdl), &s );
        int count = 0;
        MPI_Get_count( &s, mpi_datatype<decltype(recv_buf[0])>(), &count );
        recv_buf.resize(count);
      }
      MPI_Recv( recv_buf.data(), recv_buf.size(), mpi_datatype<decltype(recv_buf[0])>(), source_rank, tag, as<MPI_Comm>(hdl), &status);
      MPI_Get_count( &status, mpi_datatype<decltype(recv_buf[0])>(), &recv_count );
    } else {
      MPI_Recv( &recv_buf, 1, mpi_datatype<decltype(recv_buf)>(), source_rank, tag, as<MPI_Comm>(hdl), &status);
      MPI_Get_count( &status, mpi_datatype<decltype(recv_buf)>(), &recv_count );
    }

    return recv_count;
  }

  template <typename T>
  Request Comm::Irecv(int source_rank, T& recv_buf, int tag, bool resize_buf_with_probe ) const {
    Request req;
    MPI_Request mpi_req = MPI_REQUEST_NULL;

    if constexpr ( is_container<T>() ) {
      if ( resize_buf_with_probe ) {
        MPI_Status s;
        MPI_Probe( source_rank, tag, as<MPI_Comm>(hdl), &s );
        int count = 0;
        MPI_Get_count( &s, mpi_datatype<decltype(recv_buf[0])>(), &count );
        recv_buf.resize(count);
      }
      MPI_Irecv( recv_buf.data(), recv_buf.size(), mpi_datatype<decltype(recv_buf[0])>(), source_rank, tag, as<MPI_Comm>(hdl), &mpi_req );
    } else
      MPI_Irecv( &recv_buf, 1, mpi_datatype<decltype(recv_buf)>(), source_rank, tag, as<MPI_Comm>(hdl), &mpi_req );

    impl_cast(req.hdl).reset( &mpi_req, Raw<MPI_Request>::free );

    return req;
  }


  template < typename T >
  void Comm::broadcast( T& buffer, int root ) const {
    if constexpr ( is_container<T>() )
      MPI_Bcast( buffer.data(), buffer.size(), mpi_datatype<decltype(buffer[0])>(), root, as<MPI_Comm>(hdl) );
    else
      MPI_Bcast( &buffer, 1, mpi_datatype<decltype(buffer)>(), root, as<MPI_Comm>(hdl) );
  }


  template < typename T >
  Request Comm::Ibroadcast( T& buffer, int root ) const {
    Request req;
    MPI_Request mpi_req = MPI_REQUEST_NULL;

    if constexpr ( is_container<T>() )
      MPI_Ibcast( buffer.data(), buffer.size(), mpi_datatype<decltype(buffer[0])>(), root, as<MPI_Comm>(hdl), &mpi_req );
    else
      MPI_Ibcast( &buffer, 1, mpi_datatype<decltype(buffer)>(), root, as<MPI_Comm>(hdl), &mpi_req );

    impl_cast(req.hdl).reset( &mpi_req, Raw<MPI_Request>::free );

    return req;
  }
}
