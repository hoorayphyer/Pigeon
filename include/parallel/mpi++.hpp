#ifndef  _MPI_XX_HPP_
#define  _MPI_XX_HPP_

#include "parallel/mpi_communication.hpp"
#include "apt/pair.hpp"

namespace mpi {
  void group_free( MPI_Group* p );
  MPI_Group group_empty();

  // NOTE MPI_Group being an ordered set seems to mean that ( group_rank, process ) are stored as elements and ordered by group_rank. The mapping is solely determined by the order of processes user passed in when constructing the group.
  struct Group : apt::Handle<MPI_Group, group_free, group_empty> {
    Group() = default;
    Group ( const std::vector<int>& ranks ); // using MPI_WORLD_GROUP
    Group ( const std::vector<int>& ranks, const Group& source );

    Group& operator^= ( const Group& );
    Group& operator+= ( const Group& );
    Group& operator-= ( const Group& );

    std::vector<int> translate( const std::vector<int>& ranks, const Group& target ) const;

    int translate( int rank, const Group& target ) const;

    // translate all members in this group to the target group
    std::vector<int> translate_all( const Group& target ) const;

    int rank() const;
    int size() const;
  };

  Group operator^ ( const Group& a, const Group& b ); // intersection
  Group operator+ ( const Group& a, const Group& b ); // union
  Group operator- ( const Group& a, const Group& b ); // difference

}

namespace mpi {
  template < typename E >
  struct CommAccessor {
  private:
    inline MPI_Comm _comm() const {
      return static_cast<MPI_Comm>(static_cast<const E&>(*this));
    }
  public:
    int rank() const;
    int size() const;
    Group group() const;
  };
}

namespace mpi {
  void comm_free ( MPI_Comm* p );
  MPI_Comm comm_null();

  // a simple communicator that manages a shared memory.
  struct Comm : public apt::Handle<MPI_Comm, comm_free, comm_null>,
                public CommAccessor<Comm>,
                public P2P_Comm<Comm>,
                public Collective_Comm<Comm> {
    Comm () = default;
    // Comm ( const Comm& super_comm, const Group& sub_group ); // TODO what if sub_group is empty?



    std::optional<Comm> split ( const Group& sub_group ) const;
    std::optional<Comm> split ( std::optional<unsigned int> color, int key ) const;
  };

}

namespace mpi {
  struct CartComm : public Comm {
    CartComm( const Comm& comm, std::vector<int> dims, std::vector<bool> periodic );

    std::vector<int> rank2coords ( int rank ) const;
    int coords2rank( const std::vector<int>& coords ) const;

    inline std::vector<int> coords () const {
      return rank2coords( rank());
    }

    std::tuple<std::vector<int>, std::vector<int>>
    coords_dims() const;

    int linear_coord() const;

    apt::pair<std::optional<int>> shift(int direction, int disp = 1 ) const;
  };

}

namespace mpi {
  struct InterComm : public apt::Handle<MPI_Comm, comm_free, comm_null>,
                     public CommAccessor<Comm>,
                     public P2P_Comm<Comm>,
                     public Collective_Comm<Comm,true> {
    InterComm( const Comm& local_comm, int local_leader, const std::optional<Comm>& peer_comm, int remote_leader, int tag );

    // when peer_comm is known to all
    InterComm( const Comm& local_comm, int local_leader, const Comm& peer_comm, int remote_leader, int tag );

    int remote_size() const;
    Group remote_group() const;
  };
}

namespace mpi {
  void initialize();
  void finalize();
  extern const Comm world;
}



#endif
