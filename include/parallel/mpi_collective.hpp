#ifndef  _MPI_COLLECTIVE_HPP_
#define  _MPI_COLLECTIVE_HPP_

namespace mpi {
  enum class ReduceOp : char { SUM = 0, MAX, MAXLOC};

  template < typename Comm, typename Request >
  struct Collective_Comm {
  private:
    inline const Comm& _comm() const { return static_cast<const Comm&>(*this)._comm(); }
  public:
    void barrier() const;

    template < typename T >
    using reduce_return_t = std::optional< typename std::remove_cv<typename std::remove_reference<T>::type>::type >;

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
