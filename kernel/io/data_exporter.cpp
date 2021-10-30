#include "io/data_exporter_impl.hpp"
#include "pic.hpp"

namespace io {
using namespace pic;

template class DataExporter<real_export_t, DGrid, real_t, particle::Specs,
                            real_j_t>;
}  // namespace io
