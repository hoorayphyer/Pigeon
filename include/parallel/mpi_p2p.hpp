#ifndef  _MPI_P2P_HPP_
#define  _MPI_P2P_HPP_

namespace mpi {

  enum class SendMode : char { STD = 0, BUF, SYN, RDY };

  template < typename Comm, typename Request >
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
#endif
