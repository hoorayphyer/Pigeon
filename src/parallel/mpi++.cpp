#include "parallel/mpi++.hpp"

// global definitions
namespace mpi {
  void group_free( MPI_Group* p ) {
    if ( p && *p != MPI_GROUP_NULL )
      MPI_Group_free(p);
  }

  MPI_Group group_default() {
    return MPI_GROUP_EMPTY;
  }

  void comm_free ( MPI_Comm* p ) {
    if ( p && *p != MPI_COMM_NULL )
      MPI_Comm_free(p);
  }

  MPI_Comm comm_default() {
    return MPI_COMM_WORLD;
  }

  const Comm world;

  void initialize() {
    int is_initialized = 0;
    MPI_Initialized(&is_initialized);

    if (!is_initialized)
      MPI_Init(NULL, NULL);
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
    if ( 0 == ranks.size() ) return;
    // TODO remove duplicates
    auto* grp = new MPI_Group;
    MPI_Group_incl(source, ranks.size(), ranks.data(), grp);
    reset(grp);
  }

  Group::Group ( const std::vector<int>& ranks ) {
    if ( 0 == ranks.size() ) return;
    MPI_Group grp_world;
    MPI_Comm_group( MPI_COMM_WORLD, &grp_world);
    auto* grp = new MPI_Group;
    MPI_Group_incl(grp_world, ranks.size(), ranks.data(), grp);
    reset(grp);
    MPI_Group_free(&grp_world);
  }

  Group& Group::operator^= ( const Group& a ) {
    auto* grp = new MPI_Group;
    MPI_Group_intersection(*this, a, grp);
    reset(grp);
    return *this;
  }

  Group& Group::operator+= ( const Group& a ) {
    auto* grp = new MPI_Group;
    MPI_Group_union(*this, a, grp);
    reset(grp);
    return *this;
  }

  Group& Group::operator-= ( const Group& a ) {
    auto* grp = new MPI_Group;
    MPI_Group_difference(*this, a, grp);
    reset(grp);
    return *this;
  }

  std::vector<int> Group::translate( const std::vector<int>& ranks, const Group& target ) const {
    std::vector<int> results( ranks.size() );
    MPI_Group_translate_ranks( *this, ranks.size(), ranks.data(), target, results.data() );
    return results;
  }

  int Group::rank() const {
    int rank;
    MPI_Group_rank(*this, &rank);
    return rank;
  }

  int Group::size() const {
    int size;
    MPI_Group_size(*this, &size);
    return size;
  }

  Group operator^ ( const Group& a, const Group& b ) {
    Group result;
    auto* grp = new MPI_Group;
    MPI_Group_intersection(a, b, grp);
    result.reset(grp);
    return result;
  }

  Group operator+ ( const Group& a, const Group& b ) {
    Group result;
    auto* grp = new MPI_Group;
    MPI_Group_union(a, b, grp);
    result.reset(grp);
    return result;
  }

  Group operator- ( const Group& a, const Group& b ) {
    Group result;
    auto* grp = new MPI_Group;
    MPI_Group_difference(a, b, grp);
    result.reset(grp);
    return result;
  }

}

namespace mpi {
  template < typename E >
  int CommAccessor<E>::rank() const {
    int rank;
    MPI_Comm_rank( _comm(), &rank );
    return rank;
  };

  template < typename E >
  int CommAccessor<E>::size() const {
    int size;
    MPI_Comm_size( _comm(), &size );
    return size;
  };

  template < typename E >
  Group CommAccessor<E>::group() const {
    Group result;
    auto* grp = new MPI_Group;
    MPI_Comm_group( _comm(), grp);
    result.reset(grp);
    return result;
  }
}

// mpi::Comm
namespace mpi {
  Comm::Comm ( const Comm& super_comm, const Group& sub_group ) {
    auto* comm = new MPI_Comm;
    MPI_Comm_create_group(super_comm, sub_group, 147, comm );
    reset(comm);
  }

  Comm::Comm ( const Comm& super_comm, int color, int key ) {
    auto* comm = new MPI_Comm;
    MPI_Comm_split(super_comm, color, key, comm );
    reset(comm);
  }

}

// mpi cartesian
namespace mpi {
  CartComm::CartComm( const Comm& comm, std::vector<int> dims, std::vector<bool> periodic ) {
    const int ndims = dims.size() < periodic.size() ? dims.size() : periodic.size();

    std::vector<int> periods(ndims);
    for ( int i = 0; i < ndims; ++i )
      periods[i] = static_cast<int>( periodic[i] );

    auto* comm_cart = new MPI_Comm;
    MPI_Cart_create( comm, ndims, dims.data(), periods.data(), true, comm_cart );
    reset( comm_cart );
  }

  std::vector<int> CartComm::rank2coords ( int rank ) {
    int ndims = 0;
    MPI_Cartdim_get( *this, &ndims );
    std::vector<int> coords(ndims);
    MPI_Cart_coords( *this, rank, ndims, coords.data() );
    return coords;
  }

  int CartComm::coords2rank( const std::vector<int>& coords ) {
    int rank = 0;
    MPI_Cart_rank( *this, coords.data(), &rank);
    return rank;
  }

  int CartComm::linear_coord() {
    int result = 0;

    int ndims = 0;
    MPI_Cartdim_get( *this, &ndims );
    std::vector<int> strides(ndims+1);
    strides[0] = 1;
    std::vector<int> periods (ndims);
    std::vector<int> coords (ndims);
    MPI_Cart_get( *this, ndims, strides.data()+1, periods.data(), coords.data() );

    for ( int i = 1; i < ndims + 1; ++i )
      strides[i] *= strides[i-1];

    for ( int i = 0; i < ndims; ++i )
      result += coords[i] * strides[i];

    return result;
  }

  std::array<std::optional<int>, 2> CartComm::shift( int direction, int disp ) {
    std::array<std::optional<int>, 2> results;
    int rank_src = 0, rank_dest = 0;
    MPI_Cart_shift( *this, direction, disp, &rank_src, &rank_dest );
    if ( rank_src != MPI_PROC_NULL ) results[0].emplace(rank_src);
    if ( rank_dest != MPI_PROC_NULL ) results[1].emplace(rank_dest );
    return results;
  }

}

// intercommunicator
namespace mpi {
  InterComm::InterComm( const Comm& local_comm, int local_leader, const std::optional<Comm>& peer_comm, int remote_leader, int tag ) {
    auto* comm = new MPI_Comm;
    auto pc = ( local_comm.rank() == local_leader ) ? (*peer_comm) : MPI_COMM_NULL;
    MPI_Intercomm_create( local_comm, local_leader, pc, remote_leader, tag, comm );
    reset( comm );
  }

  int InterComm::remote_size() const {
    int size = 0;
    MPI_Comm_remote_size(*this, &size);
    return size;
  }

  Group InterComm::remote_group() const {
    Group result;

    auto* grp_remote = new MPI_Group;
    MPI_Comm_remote_group(*this, grp_remote);
    result.reset(grp_remote);

    return result;
  }

}
