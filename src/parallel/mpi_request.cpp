#include "parallel/mpi_request.hpp"

namespace mpi {
  void request_free ( MPI_Request* p ) {
    if ( p && *p != MPI_REQUEST_NULL )
      MPI_Request_free(p);
  }

  MPI_Request request_null() {
    return MPI_REQUEST_NULL;
  }

  void wait( Request& request ) {
    // MPI_Status status;
    // MPI_Wait( &request, &status );
    // return status;
    MPI_Wait( request, MPI_STATUS_IGNORE );
  }

  void waitall( std::vector<Request>& requests ) {
    MPI_Request* p = new MPI_Request [ requests.size() ];
    for ( int i = 0; i < requests.size(); ++i )
      p[i] = requests[i];
    // std::vector<MPI_Status> statuses(requests.size());
    // MPI_Waitall( requests.size(), requests.data(), statuses.data() );
    // return statuses;
    MPI_Waitall( requests.size(), p, MPI_STATUSES_IGNORE );
    delete [] p;
  }

  void cancel( Request& req ) {
    if ( static_cast<MPI_Request>(req) != MPI_REQUEST_NULL ) // this is necessary
      MPI_Cancel(req);
  }

  void cancelall( std::vector<Request>& reqs ) {
    for ( auto& req : reqs )
      cancel(req);
  }
}
