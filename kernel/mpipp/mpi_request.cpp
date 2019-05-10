#include "mpipp/mpi_request.hpp"
// RATIONALE it turns out that MPI_Wait, MPI_Waitall also manage allocation of MPI_Requests, which will interfere with apt::Handle. We want the finished mpi::Requests will be reset. In order to do that, each mpi::Request will be assigned MPI_REQUEST_NULL to avoid deallocation during reset.

// TODO MPI_Cancel still needs deallocation
// TODO persistent requests are not taken care of in wait

namespace mpi {
  void request_free ( MPI_Request* p ) {
    if ( p && *p != MPI_REQUEST_NULL )
      MPI_Request_free(p);
  }

  MPI_Request request_null() {
    return MPI_REQUEST_NULL;
  }

  void wait( Request& request ) {
    if ( !request ) return;
    MPI_Request raw = request;
    MPI_Wait( request, MPI_STATUS_IGNORE );
    if ( MPI_REQUEST_NULL == raw ) {
      MPI_Request* p = request;
      (*p) = MPI_REQUEST_NULL;
      request.reset();
    }
  }

  void waitall( std::vector<Request>& reqs ) {
    std::vector<MPI_Request> raw_reqs;
    raw_reqs.reserve(reqs.size());
    for ( int i = 0; i < reqs.size(); ++i ) {
      raw_reqs.push_back( reqs[i] ? static_cast<MPI_Request>(reqs[i]) : MPI_REQUEST_NULL );
    }

    MPI_Waitall( raw_reqs.size(), raw_reqs.data(), MPI_STATUSES_IGNORE );
    for ( int i = 0; i < reqs.size(); ++i ) {
      if ( raw_reqs[i] == MPI_REQUEST_NULL && reqs[i] ) {
        MPI_Request* p = reqs[i];
        (*p) = MPI_REQUEST_NULL;
        reqs[i].reset();
      }
    }
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
