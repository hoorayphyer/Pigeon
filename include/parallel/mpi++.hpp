#ifndef  _MPI_XX_HPP_
#define  _MPI_XX_HPP_

#include <memory>
#include <vector>
#include <optional>
#include <unordered_map>
#include <string>
#include <any>

namespace mpi {
  void initialize();
  void finalize();

  struct Handle {
  private:
    std::shared_ptr<void> _ptr {nullptr}; // type erasure
  public:
    friend std::shared_ptr<void>& impl_cast( Handle& );
    friend const std::shared_ptr<void>& impl_cast( const Handle& );
    bool operator==( Handle other ) const;
    inline bool operator!=( Handle other ) const { return !operator==(other); }
  };

  extern Handle NULL_HANDLE;

  struct Request {
    Handle hdl;
  };

  using Group = std::vector<int>;
  Group operator^ ( Group a, Group b ); // intersection
  Group operator+ ( Group a, Group b ); // union
  Group operator- ( Group a, Group b ); // difference

  enum class SendMode : char { STD = 0, BUF, SYN, RDY };
  enum class ReduceOp : char { SUM = 0, MAX, MAXLOC};

  struct Comm {
    Handle hdl;
    std::unordered_map<std::string, std::any> attrs; // this is to replace mpi raw caching mechanism

    Comm() = default;
    Comm ( Group group );

    //TODOL check assignment although it seems to work out of the box

    int rank() const;
    int size() const;

    std::vector<int> to_world( std::vector<int> ranks ) const;



    //------------------------------------
    void barrier() const;

    template <typename T_or_Container>
    void send( int dest_rank, const T_or_Container& send_buf, int tag, SendMode mode = SendMode::SYN ) const;
    template <typename T_or_Container>
    Request Isend( int dest_rank, const T_or_Container& send_buf, int tag, SendMode mode = SendMode::STD ) const;

    template <typename T_or_Container>
    int recv( int source_rank, T_or_Container& recv_buf, int tag, bool resize_buf_with_probe = false ) const; // return the actual recved number
    template <typename T_or_Container>
    Request Irecv(int source_rank, T_or_Container& recv_buf, int tag, bool resize_buf_with_probe = false ) const;


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

  extern Comm world; // TODOL make it const

  namespace topo {
    void cartesianize( Comm& comm, std::vector<int> dims, std::vector<bool> periodic );

    namespace cart {
      std::vector<int> coords ( const Comm& comm );
      std::vector<int> dims ( const Comm& comm );
      std::vector<bool> periodic ( const Comm& comm );
      int linear_coord( const Comm& comm );
    }

  }

}






#endif
