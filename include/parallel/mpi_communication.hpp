#ifndef  _MPI_COMMUNICATION_HPP_
#define  _MPI_COMMUNICATION_HPP_

#include "utility/handle.hpp"
#include <optional>

namespace mpi {
  struct Request {
    Handle hdl;
  };
}

namespace mpi {

  enum class SendMode : char { STD = 0, BUF, SYN, RDY };

  template < typename Comm >
  struct P2P_Comm {
  private:
    inline const Comm& _comm() const { return static_cast<const Comm&>(*this)._comm(); }
  public:

    template <typename T_or_Container>
    void send( int dest_rank, const T_or_Container& send_buf, int tag, SendMode mode = SendMode::SYN ) const;
    template <typename T_or_Container>
    Request Isend( int dest_rank, const T_or_Container& send_buf, int tag, SendMode mode = SendMode::STD ) const;

    template <typename T_or_Container>
    int recv( int source_rank, T_or_Container& recv_buf, int tag, bool resize_buf_with_probe = false ) const; // return the actual recved number
    template <typename T_or_Container>
    Request Irecv(int source_rank, T_or_Container& recv_buf, int tag, bool resize_buf_with_probe = false ) const;

  };

}


namespace mpi {
  enum class ReduceOp : char { SUM = 0, MAX, MAXLOC};

  template < typename T >
  using reduce_return_t = std::optional< typename std::remove_cv<typename std::remove_reference<T>::type>::type >;

  template < typename Comm >
  struct Collective_Comm {
  private:
    inline const Comm& _comm() const { return static_cast<const Comm&>(*this)._comm(); }

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
