#include "builder.hpp"

namespace pic {
template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
SimulationBuilder<DGrid, R, S, RJ, RD>::DataExporter_t& add_exporter() {
  return m_exporters.emplace_back();
}
}  // namespace pic
