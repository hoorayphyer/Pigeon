#ifndef _IO_DATA_EXPORT_HPP_
#define _IO_DATA_EXPORT_HPP_

#include <string>

struct DynamicVars;
struct Params;
namespace mpi { struct Comm; }
namespace std { template <class T> class optional; }

namespace io {
  void export_data( std::string this_run_dir,
                    int timestep,
                    const DynamicVars& dvars,
                    const Params& params,
                    const std::optional<mpi::Comm>& primary,
                    const std::optional<mpi::Comm>& ensemble );
}

#endif
