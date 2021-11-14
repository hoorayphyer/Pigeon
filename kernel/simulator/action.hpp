#pragma once
#include <string>

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

template <typename Action>
struct ActionSetter_t {
  Action& set_name(const std::string& name) {
    auto& res = static_cast<Action&>(*this);
    res.m_name = name;
    return res;
  }

  Action& set_range(std::size_t i, apt::Range range) {
    auto& res = static_cast<Action&>(*this);
    // TODO use .at(i) once we get rid of apt::array
    res.m_ranges[i] = std::move(range);
    return res;
  }
};

template <int DGrid, typename R, typename RJ>
struct FieldAction : public ActionBase<DGrid>,
                     public ActionSetter_t<FieldAction<DGrid, R, RJ>> {
  using Bundle_t = FieldBundle<DGrid, R, RJ>;
  friend class ActionSetter_t<FieldAction<DGrid, R, RJ>>;

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
}  // namespace pic
