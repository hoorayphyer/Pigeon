#include "apt/grid.hpp"
#include "apt/index.hpp"
#include "particle/shapef.hpp"
#include "simulator/argparser.hpp"
#include "simulator/builder.hpp"
#include "simulator/conf_file.hpp"

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
struct PIGEON {
  using SimulationBuilder_t = pic::SimulationBuilder<DGrid, R, S, RJ, RD>;
  using FieldAction_t = SimulationBuilder_t::FieldAction_t;
  using ParticleAction_t = SimulationBuilder_t::ParticleAction_t;
  using ConfFile_t = pic::ConfFile;
  using Grid_t = apt::Grid<R, DGrid>;
  using Index_t = apt::Index<DGrid>;
};
