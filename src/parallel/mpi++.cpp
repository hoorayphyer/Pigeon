#include "parallel/mpi++.hpp"

// global definitions
namespace mpi {
  void group_free( MPI_Group* p ) {
    if ( p && *p != MPI_GROUP_NULL )
      MPI_Group_free(p);
  }

  void comm_free ( MPI_Comm* p ) {
    if ( p && *p != MPI_COMM_NULL && *p != MPI_COMM_WORLD )
      MPI_Comm_free(p);
  }

  Comm world;

  void initialize() {
    int is_initialized = 0;
    MPI_Initialized(&is_initialized);

    if (!is_initialized)
      MPI_Init(NULL, NULL);

    // TODO check freeing deepcopied MPI_COMM_WORLD
    MPI_Comm comm_world = MPI_COMM_WORLD;
    world.reset(&comm_world);
  }

  void finalize() {
    int is_finalized = 0;
    MPI_Finalized(&is_finalized);

    if (!is_finalized) MPI_Finalize();
  }
}

// mpi::Group
namespace mpi {
  Group::Group ( const std::vector<int>& ranks, const Group& source ) {
    // TODO remove duplicates
    MPI_Group grp;
    MPI_Group_incl(source, ranks.size(), ranks.data(), &grp);
    reset(&grp);
  }

  Group::Group ( const std::vector<int>& ranks )
    : Group( ranks, world.group() ) {}

  // Group::operator std::vector<int>() {
  //   int size;
  //   MPI_Group_size(*this, &size);
  //   std::vector<int> ranks
  // }

  Group& Group::operator^= ( const Group& a ) {
    MPI_Group grp;
    MPI_Group_intersection(*this, a, &grp);
    reset(&grp);
    return *this;
  }

  Group& Group::operator+= ( const Group& a ) {
    MPI_Group grp;
    MPI_Group_union(*this, a, &grp);
    reset(&grp);
    return *this;
  }

  Group& Group::operator-= ( const Group& a ) {
    MPI_Group grp;
    MPI_Group_difference(*this, a, &grp);
    reset(&grp);
    return *this;
  }

  std::vector<int> Group::translate( const std::vector<int>& ranks, const Group& target ) const {
    std::vector<int> results( ranks.size() );
    MPI_Group_translate_ranks( *this, ranks.size(), ranks.data(), target, results.data() );
    return results;
  }

  Group operator^ ( const Group& a, const Group& b ) {
    Group result;
    MPI_Group grp;
    MPI_Group_intersection(a, b, &grp);
    result.reset(&grp);
    return result;
  }

  Group operator+ ( const Group& a, const Group& b ) {
    Group result;
    MPI_Group grp;
    MPI_Group_union(a, b, &grp);
    result.reset(&grp);
    return result;
  }

  Group operator- ( const Group& a, const Group& b ) {
    Group result;
    MPI_Group grp;
    MPI_Group_difference(a, b, &grp);
    result.reset(&grp);
    return result;
  }

}

// mpi::Comm
namespace mpi {
  Comm::Comm ( const Comm& super_comm, const Group& sub_group ) {
    MPI_Comm comm;
    MPI_Comm_create_group(super_comm, sub_group, 147, &comm );
    reset(&comm);
  }

  int Comm::rank() const {
    int rank;
    MPI_Comm_rank( *this, &rank );
    return rank;
  }

  int Comm::size() const {
    int size;
    MPI_Comm_size( *this, &size );
    return size;
  }

  Group Comm::group() const {
    Group result;
    MPI_Group grp;
    MPI_Comm_group(*this, &grp);
    result.reset(&grp);
    return result;
  }


}

// mpi cartesian
namespace mpi {
  void cartesianize( Comm& comm, std::vector<int> dims, std::vector<bool> periodic ) {
    const int ndims = dims.size() < periodic.size() ? dims.size() : periodic.size();

    std::vector<int> periods(ndims);
    for ( int i = 0; i < ndims; ++i )
      periods[i] = static_cast<int>( periodic[i] );

    MPI_Comm comm_cart;
    MPI_Cart_create( comm, ndims, dims.data(), periods.data(), true, &comm_cart );
    comm.reset( &comm_cart );
  }

  namespace cart {
    std::vector<int> rank2coords ( const Comm& comm, int rank ) {
      int ndims = 0;
      MPI_Cartdim_get( comm, &ndims );
      std::vector<int> coords(ndims);
      MPI_Cart_coords( comm, rank, ndims, coords.data() );
      return coords;
    }

    int coords2rank( const Comm& comm, const std::vector<int>& coords ) {
      int rank = 0;
      MPI_Cart_rank( comm, coords.data(), &rank);
      return rank;
    }

    int linear_coord( const Comm& comm ) {
      int result = 0;

      int ndims = 0;
      MPI_Cartdim_get( comm, &ndims );
      auto* strides = new int [ndims+1];
      strides[0] = 1;
      auto* periods = new int [ndims];
      auto* coords = new int [ndims];
      MPI_Cart_get( comm, ndims, strides+1, periods, coords );

      for ( int i = 1; i < ndims + 1; ++i )
        strides[i] *= strides[i-1];

      for ( int i = 0; i < ndims; ++i )
        result += coords[i] * strides[i];

      delete [] strides;
      delete [] periods;
      delete [] coords;

      return result;
    }

    std::array<std::optional<int>, 2> shift( const Comm& comm, int direction, int disp ) {
      std::array<std::optional<int>, 2> results;
      int rank_src = 0, rank_dest = 0;
      MPI_Cart_shift( comm, direction, disp, &rank_src, &rank_dest );
      if ( rank_src != MPI_PROC_NULL ) results[0].emplace(rank_src);
      if ( rank_dest != MPI_PROC_NULL ) results[1].emplace(rank_dest );
      return results;
    }

  }
}

// intercommunicator
namespace mpi {
  InterComm::InterComm( const Comm& local_comm, int local_leader, const std::optional<Comm>& peer_comm, int remote_leader, int tag ) {
    MPI_Comm comm;
    auto pc = ( local_comm.rank() == local_leader ) ? (*peer_comm) : MPI_COMM_NULL;
    MPI_Intercomm_create( local_comm, local_leader, pc, remote_leader, tag, &comm );
    reset( &comm );
  }

  int InterComm::remote_size() const {
    int size = 0;
    MPI_Comm_remote_size(*this, &size);
    return size;
  }

  Group InterComm::remote_group() const {
    Group result;

    MPI_Group grp_remote;
    MPI_Comm_remote_group(*this, &grp_remote);
    result.reset(&grp_remote);

    return result;
  }

}
