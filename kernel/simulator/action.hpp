#pragma once
#include <string>

#include "apt/range.hpp"
#include "simulator/bundle.hpp"

namespace pic {
template <int DGrid>
struct ActionBase : public array<apt::Range, DGrid> {
 private:
  std::string m_name = "Unknown";

 public:
  void setName(std::string name) { m_name = std::move(name); }

  const auto& name() const noexcept { return m_name; }

  virtual ~ActionBase() = default;
};

template <typename R, int DGrid, typename RJ>
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
}  // namespace pic
