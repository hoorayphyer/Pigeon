#include "apt/grid.hpp"
#include "apt/index.hpp"
#include "particle/array.hpp"
#include "particle/shapef.hpp"
#include "particle/species_predef.hpp"
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

  template <typename ConcreteAction>
  using InitialConditionAction_t = pic::ActionWithSetters<
      ConcreteAction, DGrid,
      typename SimulationBuilder_t::InitialConditionAction_t>;

  template <typename ConcreteAction>
  using PostResumeAction_t =
      pic::ActionWithSetters<ConcreteAction, DGrid,
                             typename SimulationBuilder_t::PostResumeAction_t>;

  using ExportBundle_t = SimulationBuilder_t::ExportBundle_t;

  using ConfFile_t = pic::ConfFile;
  using Grid_t = apt::Grid<R, DGrid>;
  using Index_t = apt::Index<DGrid>;
  using Vec3 = apt::Vec<R, 3>;
  using Particle = particle::Particle<R, S>;
  using vParticle = particle::vParticle<R, S>;

  using ParticleArray_t = particle::array<R, S>;

  template <int DField>
  using Field = field::Field<R, DField, DGrid>;

  using JField = field::Field<RJ, 3, DGrid>;  // TODO get rid of this

  using IOField = field::Field<RD, 3, DGrid>;
  using IOGrid = apt::Grid<RD, DGrid>;  // TODO can we just use same grid???

  using Schedule = pic::Schedule;
  using ExportSchedule = pic::ExportSchedule;
  using CheckpointSchedule = pic::CheckpointSchedule;
  using DynamicLoadBalanceSchedule = pic::DynamicLoadBalanceSchedule;
  using ProfilingSchedule = pic::ProfilingSchedule;
};

using particle::species;
