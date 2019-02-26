#ifndef  _MPI_COMMUNICATION_HPP_
#define  _MPI_COMMUNICATION_HPP_

#include "apt/handle.hpp"
#include <mpi.h>
#include <type_traits>
#include <optional>
#include <vector>

namespace mpi {
  template <typename T>
  constexpr MPI_Datatype datatype(const T* = nullptr) noexcept {
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

  template <typename T>
  constexpr MPI_Datatype datatype(const T&) noexcept {
    return datatype((const T*)0);
  }

}


namespace mpi {
  void request_free ( MPI_Request* p );
  MPI_Request request_default();
  using Request = apt::Handle<MPI_Request, request_free, request_default >;
}

namespace std {
  template < class, class > class vector;
}

namespace mpi {

  enum class SendMode : char { STD = 0, BUF, SYN, RDY };

  void wait( Request& req );
  void waitall( std::vector<Request>& reqs );

  void cancel( Request& req );
  void cancelall( std::vector<Request>& reqs );


  template < typename Comm >
  struct P2P_Comm {
  private:
    inline MPI_Comm _comm() const {
      return static_cast<MPI_Comm>(static_cast<const Comm&>(*this));
    }

  public:
    template < typename T >
    int probe( int source_rank, int tag, const T* ) const; // returns the send count

    template < typename T >
    inline int probe( int source_rank, int tag, const T& ) const {
      return probe( source_rank, tag, (const T*)0 );
    }

    template <typename T>
    void send( int dest_rank, int tag, const T* send_buf, int send_count = 1, SendMode mode = SendMode::SYN ) const;
    template <typename T>
    Request Isend( int dest_rank, int tag, const T* send_buf, int send_count = 1, SendMode mode = SendMode::STD ) const;

    template <typename T>
    int recv( int source_rank, int tag , T* recv_buf, int recv_count_max) const; // return the actual recved number
    template <typename T>
    Request Irecv(int source_rank, int tag, T* recv_buf, int recv_count_max ) const;

  };

}


namespace mpi {
  enum class ReduceOp : char { SUM = 0, MAX, MAXLOC};

  template < typename T >
  using reduce_return_t = std::optional< typename std::remove_cv<typename std::remove_reference<T>::type>::type >;

  template < typename Comm >
  struct Collective_Comm {
  private:
    inline MPI_Comm _comm() const {
      return static_cast<MPI_Comm>(static_cast<const Comm&>(*this));
    }

  public:
    void barrier() const;

    template < bool In_Place, ReduceOp Op, typename T_or_Container>
    reduce_return_t<T_or_Container> reduce( T_or_Container& buffer, int root ) const;
    template < bool In_Place, ReduceOp Op, typename T_or_Container>
    std::tuple< Request, reduce_return_t<T_or_Container> > Ireduce( T_or_Container& buffer, int root ) const;

    template < typename T_or_Container >
    void broadcast( T_or_Container& buffer, int root ) const;
    template < typename T_or_Container >
    Request Ibroadcast( T_or_Container& buffer, int root ) const;


    // // recvcount refers to counts of recved data from one single process, rather
    // // than the total counts from all processes.
    // // supports passing in MPI_IN_PLACE as the send_buf at root
    // // use RecvTypePtr to allow use of nullptr at non-root processes
    // template <typename SendType, typename RecvTypePtr>
    // void gather( const SendType* send_buf, int sendcount, RecvTypePtr recv_buf, int recvcount, int root ) const;

    // template <typename SendType, typename RecvTypePtr>
    // MPI_Request Igather( const SendType* send_buf, int sendcount, RecvTypePtr recv_buf, int recvcount, int root ) const;

    // // variable gather
    // template <typename SendType, typename RecvTypePtr>
    // void gatherv( const SendType* send_buf, int sendcount, RecvTypePtr recv_buf, const int* recvcounts, const int* displs, int root ) const;

    // template <typename SendType, typename RecvTypePtr>
    // MPI_Request Igatherv( const SendType* send_buf, int sendcount, RecvTypePtr recv_buf, const int* recvcounts, const int* displs, int root ) const;

  };

}
#endif
