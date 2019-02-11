#ifndef  _MPI_XX_HPP_
#define  _MPI_XX_HPP_

#include "parallel/mpi_communication.hpp"
#include <vector>
#include <array>

namespace mpi {
  void initialize();
  void finalize();

  using Group = std::vector<int>;
  Group operator^ ( Group a, Group b ); // intersection
  Group operator+ ( Group a, Group b ); // union
  Group operator- ( Group a, Group b ); // difference

  struct Comm : public P2P_Comm<Comm>,
                public Collective_Comm<Comm>{
  private:
    const Comm& _comm() const noexcept { return *this; }
  public:
    friend class P2P_Comm<Comm>;
    friend class Collective_Comm<Comm>;
    Handle hdl;

    Comm() = default;
    Comm ( Group group );

    //TODOL check assignment although it seems to work out of the box

    int rank() const;
    int size() const;

    std::vector<int> to_world( std::vector<int> ranks ) const;

  };

  extern Comm world; // TODOL make it const
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
  struct InterComm : public P2P_Comm<InterComm>,
                     public Collective_Comm<InterComm> {
  private:
    const InterComm& _comm() const noexcept { return *this; }
  public:
    friend class P2P_Comm<InterComm>;
    friend class Collective_Comm<InterComm>;

    Handle hdl;

    InterComm( const Comm& local_comm, int local_leader, const std::optional<Comm>& peer_comm, int remote_leader, int tag );
    int remote_size() const;
    Group remote_group() const;
  };
}





#endif
