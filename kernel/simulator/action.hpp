#pragma once
#include <functional>
#include <string>
#include <type_traits>

#include "apt/range.hpp"
#include "simulator/bundle.hpp"

namespace pic {
template <int DGrid>
struct ActionBase {
 protected:
  apt::array<apt::Range, DGrid> m_ranges;
  std::string m_name = "Unknown";

 public:
  const auto& name() const noexcept { return m_name; }

  const auto& ranges() const noexcept { return m_ranges; }

  auto& ranges() noexcept { return m_ranges; }

  virtual ~ActionBase() = default;
};

template <int DGrid, typename R, typename RJ>
struct FieldAction : public ActionBase<DGrid> {
  using Bundle_t = FieldBundle<DGrid, R, RJ>;

  virtual ~FieldAction() = default;

  virtual void operator()(const Bundle_t& bundle) const = 0;
};

template <int DGrid, typename R, template <typename> class S, typename RJ>
struct ParticleAction : public ActionBase<DGrid> {
  using Bundle_t = ParticleBundle<DGrid, R, S, RJ>;
  virtual ~ParticleAction() = default;

  virtual void operator()(const Bundle_t& bundle) const = 0;
};

template <int DGrid, typename R, template <typename> class S, typename RJ>
struct InitialConditionAction : public ActionBase<DGrid> {
  using Bundle_t = InitialConditionBundle<DGrid, R, S, RJ>;
  virtual ~InitialConditionAction() = default;

  virtual void operator()(const Bundle_t& bundle) const = 0;
};

template <int DGrid, typename R, template <typename> class S, typename RJ>
struct PostResumeAction : public ActionBase<DGrid> {
  using Bundle_t = PostResumeBundle<DGrid, R, S, RJ>;
  virtual ~PostResumeAction() = default;

  virtual void operator()(const Bundle_t& bundle) const = 0;
};

// RATIONALE
//
// want to add a few setters in a ConcreteAction class and return
// ConcreteAction&, so as to be chained with other setters that class may
// have.
//
// AbstractAction is one of the above Action types that derive from
// ActionBase. The setters only set members of ActionBase, but adding
// AbstractAction and automatically inheriting AbstractAction makes the user
// code simpler. For example,
//
//     struct A : ActionWithSetter<A, DGrid, FieldAction>;
//
// `A` automatically inherits FieldAction while getting the setters.
// Otherwise, `A` has to add a second inheritance to FieldAction to register
// itself as a FieldAction.
//
// Why not just make FieldAction template on ConcreteAction and put these
// setters therein? Well, we need a common base class for all concrete field
// actions, so as to be put in the SimulationBuilder.
template <typename ConcreteAction, int DGrid, typename AbstractAction>
struct ActionWithSetters : public AbstractAction {
  static_assert(std::is_base_of_v<ActionBase<DGrid>, AbstractAction>);

  ConcreteAction& set_name(const std::string& name) {
    ActionBase<DGrid>::m_name = name;
    return static_cast<ConcreteAction&>(*this);
  }

  ConcreteAction& set_range(std::size_t i, apt::Range range) {
    // TODO use .at(i) once we get rid of apt::array
    ActionBase<DGrid>::m_ranges[i] = std::move(range);
    return static_cast<ConcreteAction&>(*this);
  }

  ConcreteAction& apply(std::function<void(ConcreteAction&)> f) {
    auto& act = static_cast<ConcreteAction&>(*this);
    f(act);
    return act;
  }
};
}  // namespace pic
