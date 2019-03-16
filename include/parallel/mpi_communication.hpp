#ifndef  _MPI_COMMUNICATION_HPP_
#define  _MPI_COMMUNICATION_HPP_

#include "parallel/mpi_datatype.hpp"
#include "apt/handle.hpp"
#include <mpi.h>
#include <optional>
#include <vector>

namespace mpi {
  void request_free ( MPI_Request* p );
  MPI_Request request_null();
  using Request = apt::Handle<MPI_Request, request_free, request_null >;
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
  enum class by : char { SUM = 0, MAX, MAXLOC};

  constexpr bool INPLACE = true;

  template < typename Comm, bool Inter = false >
  struct Collective_Comm {
  private:
    inline MPI_Comm _comm() const {
      return static_cast<MPI_Comm>(static_cast<const Comm&>(*this));
    }

  public:
    void barrier() const;

    template < by Op, bool In_Place = false, typename T >
    std::optional<std::vector<T>>
    reduce( T* buffer, int count, int root ) const;

    // TODO fix this. Returning optional on nonblocking call may not make sense
    // template < by Op, bool In_Place = false, typename T >
    // std::tuple< Request, std::optional<std::vector<T> > >
    // Ireduce( T* buffer, int count, int root ) const;

    template < typename T >
    void broadcast( int root, T* buffer, int count ) const;

    template < typename T >
    Request Ibroadcast( int root, T* buffer, int count ) const;

    template < typename T >
    std::vector<T> allgather( const T* send_buf, int send_count ) const;

    // scatter for intra_comm is assumed to be in_place. This means in any case, root will not send anything to itself
    template < typename T >
    void scatter( int root, T* buffer, int count ) const;
  };

}
#endif
