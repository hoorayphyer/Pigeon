#pragma once

#include <mpi.h>

#include <vector>

#include "apt/handle.hpp"

namespace mpi {
void request_free(MPI_Request* p);
MPI_Request request_null();
using Request = apt::Handle<MPI_Request, request_free, request_null>;

// NOTE requests will be reset when finished waiting
void wait(Request& req);
void waitall(std::vector<Request>& reqs);

void cancel(Request& req);
void cancelall(std::vector<Request>& reqs);
}  // namespace mpi
