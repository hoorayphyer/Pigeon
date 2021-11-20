#include "apt/grid.hpp"
#include "apt/index.hpp"
#include "particle/array.hpp"
#include "particle/shapef.hpp"
#include "particle/species_predef.hpp"
#include "simulator/action_predefined.hpp"
#include "simulator/argparser.hpp"
#include "simulator/builder.hpp"
#include "simulator/conf_file.hpp"

template <int DGrid, typename R, template <typename> class S, typename RJ,
          typename RD>
struct PIGEON {
  using SimulationBuilder_t = pic::SimulationBuilder<DGrid, R, S, RJ, RD>;

  template <typename ConcreteAction>
  using FieldAction_t =
      pic::ActionWithSetters<ConcreteAction, DGrid,
                             typename SimulationBuilder_t::FieldAction_t>;

  template <typename ConcreteAction>
  using ParticleAction_t =
      pic::ActionWithSetters<ConcreteAction, DGrid,
                             typename SimulationBuilder_t::ParticleAction_t>;

  using ParticleMigrator_t = particle::Migrator<DGrid, R, S, RJ>;

  using ConfFile_t = pic::ConfFile;
  using Grid_t = apt::Grid<R, DGrid>;
  using Index_t = apt::Index<DGrid>;
  using Vec3 = apt::Vec<R, 3>;
  using Particle = particle::Particle<R, S>;

  using ParticleArray_t = particle::array<R, S>;

  template <int DField>
  using Field = field::Field<R, DField, DGrid>;
};

using particle::species;
