#ifndef _MPI_POINT_TO_POINT_HPP_
#define _MPI_POINT_TO_POINT_HPP_

#include "parallel/mpi_request.hpp"

namespace mpi {

  enum class SendMode : char { STD = 0, BUF, SYN, RDY };

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

#endif
