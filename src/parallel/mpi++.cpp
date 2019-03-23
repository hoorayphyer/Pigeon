#include "parallel/mpi++.hpp"

// global definitions
namespace mpi {
  void group_free( MPI_Group* p ) {
    if ( p && *p != MPI_GROUP_NULL )
      MPI_Group_free(p);
  }

  MPI_Group group_empty() {
    return MPI_GROUP_EMPTY;
  }

  void comm_free ( MPI_Comm* p ) {
    if ( p && *p != MPI_COMM_NULL )
      MPI_Comm_free(p);
  }

  MPI_Comm comm_null() {
    return MPI_COMM_NULL;
  }

  Comm make_world() {
    Comm comm;
    comm.set_fallback_handle(MPI_COMM_WORLD);
    return comm;
  }
  const Comm world = make_world();

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

  int Group::translate( int rank, const Group& target ) const {
    int res;
    MPI_Group_translate_ranks( *this, 1, &rank, target, &res );
    return res;
  }

  std::vector<int> Group::translate_all( const Group& target ) const {
    std::vector<int> ranks( size() );
    for ( int i = 0; i < size(); ++i ) ranks[i] = i;
    return translate( ranks, target );
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
  // Comm::Comm ( const Comm& super_comm, const Group& sub_group ) {
  //   auto* comm = new MPI_Comm;
  //   MPI_Comm_create_group(super_comm, sub_group, 147, comm );
  //   reset(comm);
  // }

  std::optional<Comm> Comm::split ( const Group& sub_group ) const {
    std::optional<Comm> res;
    auto* comm = new MPI_Comm;
    MPI_Comm_create(*this, sub_group, comm );
    if ( *comm != MPI_COMM_NULL ) {
      res.emplace();
      res->reset(comm);
    } else {
      delete comm;
    }
    return res;
  }

  std::optional<Comm> Comm::split ( std::optional<unsigned int> color, int key ) const {
    std::optional<Comm> res;
    auto* comm = new MPI_Comm;
    MPI_Comm_split(*this, color ? *color : MPI_UNDEFINED, key, comm );
    if ( *comm != MPI_COMM_NULL ) {
      res.emplace();
      res->reset(comm);
    } else {
      delete comm;
    }
    return res;
  }
}

// mpi cartesian
namespace mpi {
  CartComm::CartComm( const Comm& comm, std::vector<int> dims, std::vector<bool> periodic ) {
    // TODOL enforce same dimensionality through interface, for example, use std::vector<apt::pair>
    const int ndims = dims.size() < periodic.size() ? dims.size() : periodic.size();

    std::vector<int> periods(ndims);
    for ( int i = 0; i < ndims; ++i )
      periods[i] = static_cast<int>( periodic[i] );

    auto* comm_cart = new MPI_Comm;
    MPI_Cart_create( comm, ndims, dims.data(), periods.data(), true, comm_cart );
    reset( comm_cart );
  }

  std::vector<int> CartComm::rank2coords ( int rank ) const {
    int ndims = 0;
    MPI_Cartdim_get( *this, &ndims );
    std::vector<int> coords(ndims);
    MPI_Cart_coords( *this, rank, ndims, coords.data() );
    return coords;
  }

  int CartComm::coords2rank( const std::vector<int>& coords ) const {
    int rank = 0;
    MPI_Cart_rank( *this, coords.data(), &rank);
    return rank;
  }

  std::tuple<std::vector<int>, std::vector<int>, std::vector<bool> >
  CartComm::coords_dims_periodic() const {
    int ndims = 0;
    MPI_Cartdim_get( *this, &ndims );
    std::vector<int> coords (ndims);
    std::vector<int> dims(ndims);
    std::vector<int> is_per (ndims);
    MPI_Cart_get( *this, ndims, dims.data(), is_per.data(), coords.data() );
    std::vector<bool> periodic(ndims);
    for ( int i = 0; i < ndims; ++i )
      periodic[i] = is_per[i];
    return make_tuple( coords, dims, periodic );
  }

  apt::pair<std::optional<int>> CartComm::shift( int ith_dim, int disp ) const {
    apt::pair<std::optional<int>> results;
    int rank_src = 0, rank_dest = 0;
    MPI_Cart_shift( *this, ith_dim, disp, &rank_src, &rank_dest );
    if ( rank_src != MPI_PROC_NULL ) results[LFT].emplace(rank_src);
    if ( rank_dest != MPI_PROC_NULL ) results[RGT].emplace(rank_dest );
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

  InterComm::InterComm( const Comm& local_comm, int local_leader, const mpi::Comm& peer_comm, int remote_leader, int tag ) {
    auto* comm = new MPI_Comm;
    MPI_Intercomm_create( local_comm, local_leader, peer_comm, remote_leader, tag, comm );
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
