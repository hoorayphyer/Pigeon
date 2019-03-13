#ifndef  _DOMAIN_H_
#define  _DOMAIN_H_

#include "Fields.h"
#include <mpi.h>
#include "ArrayOperations.h"
#include "FUParams.h"

enum CommTags {
  SEND_LEFT = 0,
  SEND_RIGHT
};

struct MPICartComm {
private:
  const MPI_Comm _comm;

  void Handle_MPI_Error(int error_code) const {
    constexpr int BUFSIZE = 1024;
    if (error_code != MPI_SUCCESS) {
      char error_string[BUFSIZE];
      int length_of_error_string;

      int rank = 0;
      MPI_Comm_rank( _comm, &rank );

      MPI_Error_string(error_code, error_string, &length_of_error_string);
      fprintf(stderr, "%s rank %3d: %s\n", "old_field_solver", rank, error_string);
    }
  }

  template < typename T >
  constexpr MPI_Datatype get_mpi_datatype ( T ) const noexcept {
    if constexpr ( std::is_same_v<T, char> ) return MPI_CHAR;
    else if ( std::is_same_v<T, int> ) return MPI_INT;
    else if ( std::is_same_v<T, float> ) return MPI_FLOAT;
    else if ( std::is_same_v<T, double> ) return MPI_DOUBLE;
    else return MPI_DATATYPE_NULL;
  }

public:
  MPICartComm( const MPI_Comm& comm ) : _comm(comm) {}

  template <typename T>
  MPI_Request Isend(int dest_rank, int tag, const T* values, int n ) const {

    MPI_Request request;
    MPI_Datatype type = get_mpi_datatype(*values);

    int error_code = MPI_Isend((void*)values, n, type, dest_rank, tag, _comm, &request);

    Handle_MPI_Error(error_code);
    return request;
  }

  template <typename T>
  MPI_Request Irecv(int source_rank, int tag, T* values, int n) const {
    MPI_Request request;
    MPI_Datatype type = get_mpi_datatype(*values);

    int error_code = MPI_Irecv((void*)values, n, type, source_rank, tag, _comm, &request);

    Handle_MPI_Error(error_code);
    return request;
  }

};


struct MPIReqPair{
private:
  MPI_Request _reqs[2] { MPI_REQUEST_NULL, MPI_REQUEST_NULL };
public:
  ~MPIReqPair() {
    for ( int i = 0; i < 2; ++i ) {
      if (_reqs[i] != MPI_REQUEST_NULL )
        MPI_Request_free(_reqs+i);
    }
  }

  inline MPI_Request& operator[]( bool i ) noexcept { return _reqs[i]; }

  inline void waitall() {
    MPI_Waitall( 2, _reqs, MPI_STATUSES_IGNORE );
  }
};

class FieldCommunicator {
  typedef MultiArray<Scalar> array_type;

  const MPICartComm& _comm;
  const FUParams& _params;

  std::array<array_type, 3> _bufferFieldRecv; ///< Receive buffer for field guard cell exchange
  std::array<array_type, 3> _bufferFieldSend; ///< Send buffer for field guard cell exchange

  template <typename T>
  void SendGuardCellsLeftRight(int direction, CommTags leftright, MultiArray<T>& array, const Grid& grid) {
    if (direction >= 3 || direction < 0)
      throw std::invalid_argument("Invalid direction!");

    // Obtain the starting index of send and receive buffers in the grid
    Index sendId, recvId;
    Extent sendExt;
    if (direction == 0) {
      sendId = Index(
          (leftright == SEND_LEFT ? grid.guard[0] : grid.reducedDim(0)), 0, 0);
      recvId = Index(
          (leftright == SEND_LEFT ? grid.dims[0] - grid.guard[0] : 0), 0, 0);
      sendExt = Extent(grid.guard[0], grid.dims[1], grid.dims[2]);
    } else if (direction == 1) {
      sendId = Index(
          0, (leftright == SEND_LEFT ? grid.guard[1] : grid.reducedDim(1)), 0);
      recvId = Index(
          0, (leftright == SEND_LEFT ? grid.dims[1] - grid.guard[1] : 0), 0);
      sendExt = Extent(grid.dims[0], grid.guard[1], grid.dims[2]);
    } else if (direction == 2) {
      sendId = Index(
          0, 0, (leftright == SEND_LEFT ? grid.guard[2] : grid.reducedDim(2)));
      recvId = Index(
          0, 0, (leftright == SEND_LEFT ? grid.dims[2] - grid.guard[2] : 0));
      sendExt = Extent(grid.dims[0], grid.dims[1], grid.guard[2]);
    }

    const auto& neighbor_left = _params.neighbor_left;
    const auto& neighbor_right = _params.neighbor_right;
    // Determine the from and destination rank
    int rank_from = (leftright == SEND_LEFT ? neighbor_right[direction] : neighbor_left[direction]);
    int rank_dest = (leftright == SEND_LEFT ? neighbor_left[direction] : neighbor_right[direction]);

    MPIReqPair reqs;

    if ( rank_dest != NEIGHBOR_NULL ) {
      // Copy the content to the send buffer
      copy_to_linear( _bufferFieldSend[direction].begin(), array.index(sendId), sendExt );
      reqs[0] = _comm.Isend(rank_dest, leftright, _bufferFieldSend[direction].data(), sendExt.size());
    }

    if ( rank_from != NEIGHBOR_NULL ) {
      reqs[1] = _comm.Irecv(rank_from, leftright, _bufferFieldRecv[direction].data(), sendExt.size());
    }

    // wait before add_from_linear
    reqs.waitall();

    if ( rank_from != NEIGHBOR_NULL ) {
      copy_from_linear(array.index(recvId), _bufferFieldRecv[direction].begin(), sendExt);
    }

  }

  template <typename T>
  void SendGuardCellsLeftRight(int direction, CommTags leftright, MultiArray<T>& array) {
    SendGuardCellsLeftRight(direction, leftright, array, _params.grid);
  }
  // void SendGuardCellsRight(int direction, MultiArray<Scalar, DIM>& array);


  // copy add information in guard cells to bulk. Used in current deposition
  template <typename T>
  void SendAddCellsLeftRight(int direction, CommTags leftright, MultiArray<T>& array,
                             const Index& sendId, const Index& recvId, const Extent& sendExt) {
    if (direction >= 3 || direction < 0)
      throw std::invalid_argument("Invalid direction!");

    const auto& neighbor_left = _params.neighbor_left;
    const auto& neighbor_right = _params.neighbor_right;
    // Determine the from and destination rank
    int rank_from = ( leftright == SEND_LEFT ? neighbor_right[direction] : neighbor_left[direction] );
    int rank_dest = ( leftright == SEND_LEFT ? neighbor_left[direction] : neighbor_right[direction] );

    MPIReqPair reqs;
    if ( rank_dest != NEIGHBOR_NULL ) {
      // Copy the content to the send buffer
      copy_to_linear( _bufferFieldSend[direction].begin(), array.index(sendId), sendExt );
      // Erase sent-away contents
      fill( array.index(sendId), sendExt, 0.0 );
      reqs[0] = _comm.Isend(rank_dest, leftright, _bufferFieldSend[direction].data(), sendExt.size());
    }

    if ( rank_from != NEIGHBOR_NULL ) {
      reqs[1] = _comm.Irecv(rank_from, leftright, _bufferFieldRecv[direction].data(), sendExt.size());
    }

    // wait before add_from_linear
    reqs.waitall();

    if ( rank_from != NEIGHBOR_NULL ) {
      add_from_linear(array.index(recvId), _bufferFieldRecv[direction].begin(), sendExt);
    }

  }

public:
  FieldCommunicator(const MPI_Comm& comm, const FUParams& params )
    : _comm(comm), _params(params) {

    const auto& grid = params.grid;

    // Initialize field send and recv buffers
    _bufferFieldRecv[0] = array_type(grid.guard[0], grid.dims[1], grid.dims[2]);
    _bufferFieldSend[0] = array_type(grid.guard[0], grid.dims[1], grid.dims[2]);
    _bufferFieldRecv[0].assign(0.0);
    _bufferFieldSend[0].assign(0.0);
    if (grid.dimension >= 2) {
      _bufferFieldRecv[1] = array_type(grid.dims[0], grid.guard[1], grid.dims[2]);
      _bufferFieldSend[1] = array_type(grid.dims[0], grid.guard[1], grid.dims[2]);
      _bufferFieldRecv[1].assign(0.0);
      _bufferFieldSend[1].assign(0.0);
    }
    if (grid.dimension >= 3) {
      _bufferFieldRecv[2] = array_type(grid.dims[0], grid.dims[1], grid.guard[2]);
      _bufferFieldSend[2] = array_type(grid.dims[0], grid.dims[1], grid.guard[2]);
      _bufferFieldRecv[2].assign(0.0);
      _bufferFieldSend[2].assign(0.0);
    }

  }

  // overloads of SendGuardCells
  template <typename T>
  void SendGuardCells(ScalarField<T>& field) {
    SendGuardCells(field.data(), field.grid());
  }

  template <typename T>
  void SendGuardCells(VectorField<T>& field, int dir) {
    for (int i = 0; i < VECTOR_DIM; ++i) {
      SendGuardCells(field.data(i), field.grid(), dir);
    }
  }

  template <typename T>
  void SendGuardCells(VectorField<T>& field) {
    for (int i = 0; i < VECTOR_DIM; ++i) {
      SendGuardCells(field.data(i), field.grid());
    }
  }

  template <typename T>
  void SendGuardCells(MultiArray<T>& field, const Grid& grid, int dir) {
    if ( dir < grid.dimension ) {
      // First send left, then send right, in direction dir
      SendGuardCellsLeftRight(dir, SEND_LEFT, field, grid);
      SendGuardCellsLeftRight(dir, SEND_RIGHT, field, grid);
    }
  }

  template <typename T>
  void SendGuardCells(MultiArray<T>& field, const Grid& grid) {
    for (int i = 0; i < grid.dimension; ++i) {
      // First send left, then send right, in direction i
      SendGuardCellsLeftRight(i, SEND_LEFT, field, grid);
      SendGuardCellsLeftRight(i, SEND_RIGHT, field, grid);
    }
  }

  template <typename T>
  void SendGuardCells(MultiArray<T>& field) {
    SendGuardCells(field, _params.grid);
  }

  // SendAdd in longitudinal and transverse directions
  template < typename T >
  void SendAddCells ( MultiArray<T>& array, const Grid& grid, int dir, bool scanLong = false ) {
    if ( dir >= grid.dimension ) return;

    Index sendId( 0, 0, 0 );
    Index recvId( 0, 0, 0 );
    Index sendExt( grid.dims[0], grid.dims[1], grid.dims[2] );
    sendExt[dir] = grid.guard[dir];

    //first send right
    //Note: it is crucial that sentaway values are cleared in the case of scanning in the longitudinal.
    //FIXME: the case of unstaggered array is not treated yet.
    sendId[dir] = grid.dims[dir] - grid.guard[dir] - static_cast<int>(scanLong);
    recvId[dir] = grid.guard[dir] - static_cast<int>(scanLong);
    SendAddCellsLeftRight( dir, SEND_RIGHT, array, sendId, recvId, sendExt );

    //then send left
    sendId[dir] = 0;
    recvId[dir] = grid.reducedDim(dir);
    SendAddCellsLeftRight( dir, SEND_LEFT, array, sendId, recvId, sendExt );

    return;
  }

};


#endif   // ----- #ifndef _DOMAIN_H_  -----


