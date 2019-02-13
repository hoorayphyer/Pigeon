#ifndef  _MPI_XX_HPP_
#define  _MPI_XX_HPP_

#include "parallel/mpi_communication.hpp"
#include <array>

namespace mpi {
  void group_free( MPI_Group* p ) {
    if ( p && *p != MPI_GROUP_NULL )
      MPI_Group_free(p);
  }

  // TODO so far, no use of Group. make group look like an ordered set
  struct Group : apt::Handle<MPI_Group, group_free> {
    // Group operator^ ( const Group& a, const Group& b ); // intersection
    // Group operator+ ( const Group& a, const Group& b ); // union
    // Group operator- ( const Group& a, const Group& b ); // difference
  };

}

namespace mpi {
  void comm_free ( MPI_Comm* p ) {
    if ( p && *p != MPI_COMM_NULL && *p != MPI_COMM_WORLD )
      MPI_Comm_free(p);
  }

  struct Comm : public apt::Handle<MPI_Comm, comm_free>,
                public P2P_Comm<Comm>,
                public Collective_Comm<Comm> {
  private:
    const Comm& _comm() const noexcept { return *this; }
  public:
    friend class P2P_Comm<Comm>;
    friend class Collective_Comm<Comm>;

    Comm() = default;
    Comm( MPI_Comm* mpi_comm );
    Comm ( const Group& group );

    int rank() const;
    int size() const;

    std::vector<int> to_world( std::vector<int> ranks ) const;

  };

}

namespace mpi {
  void cartesianize( Comm& comm, std::vector<int> dims, std::vector<bool> periodic );

  namespace cart {
    std::vector<int> rank2coords ( const Comm& comm, int rank );
    int coords2rank( const Comm& comm, const std::vector<int>& coords );

    inline std::vector<int> coords ( const Comm& comm ) {
      return rank2coords(comm, comm.rank());
    }

    int linear_coord( const Comm& comm );

    std::array<std::optional<int>, 2> shift( const Comm& comm, int direction, int disp );

  }
}

namespace mpi {
  struct InterComm : public Comm {
    InterComm( const Comm& local_comm, int local_leader, const std::optional<Comm>& peer_comm, int remote_leader, int tag );
    int remote_size() const;
    Group remote_group() const;
  };
}

namespace mpi {
  void initialize();
  void finalize();
  extern Comm world; // TODOL make it const
}



#endif
