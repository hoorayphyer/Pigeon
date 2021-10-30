#ifndef _PIC_TRACING_HPP_
#define _PIC_TRACING_HPP_

#include <unordered_map>

#include "apt/bit_manip.hpp"
#include "particle/action.hpp"
#include "particle/array.hpp"
#include "particle/map.hpp"
#include "particle/particle.hpp"
#include "particle/properties.hpp"
#include "pic/plans.hpp"

// Tracing Logic
// 1. a particle is uniquely identified by (sp, wr, sn), which we call id. wr is
// used in favor of ens_label because wr stays same throughout
// 2. the id is assigned only when a particle is first time traced.
// 3. sn begins with 1. sn == 0 is reversed to signal that the particle has
// never been traced before.
// 4. there is a flag::traced dedicated to tracing. When tracing, this flag is
// set, and id must be assigned accordingly. When untracing, simply reset this
// flag, leaving id untouched so that later another tracing doesn't assign new
// id to this particle
// 5. it is considered inconsistent if a particle has flag:traced set yet has sn
// == 0. Implementation should prevent this. In particular, setting flag::traced
// through state should be disabled
// 6. NOTE one problem might be the serial number, once assigned, cannot by
// recycled. So if in the future, even if user decides to ignore part of the
// traced set, the counter is determined by the max count remaining. NOTE adding
// yet another dimension, say batch No., may mitigate this issue, but it has its
// own limitations, such as limited max number of particles that can be traced
// within one batch, assuming same number of bits used.

namespace particle {
struct TracingManager {
 private:
  static constexpr int Nbits_wr = 12;
  static constexpr int Nbits_sn = 23;
  static_assert(Nbits_wr > 0 and Nbits_sn > 0);

 public:
  static void init(std::size_t wr,
                   std::unordered_map<species, std::size_t> counter) {
    _data()._wr = wr;
    _data()._counter = std::move(counter);
  }

  template <typename R, template <typename> class S>
  static void init(std::size_t wr, const map<array<R, S>>& particles) {
    return init(wr, parse_tracing(particles));
  }

  template <typename P>
  static void trace(P&& ptc) {  // TODOL sematics, preferably P&
    std::size_t id = ptc.template get<pid>();
    if (apt::getbits<0, 1, bool>(id)) return;

    apt::setbits<0, 1>(id, true);
    auto sn = apt::getbits<1 + Nbits_wr, Nbits_sn>(id);
    if (sn == 0) {
      sn = ++_data()._counter[ptc.template get<species>()];
      apt::setbits<1, Nbits_wr>(id, _data()._wr);
      apt::setbits<1 + Nbits_wr, Nbits_sn>(id, sn);
    }
    ptc.template set<pid>(id);
  }

  template <typename P>
  static void untrace(P&& ptc) {  // TODOL sematics, preferably P&
    std::size_t id = ptc.template get<pid>();
    apt::setbits<0, 1>(id, false);
    ptc.template set<pid>(id);
  }

  template <typename P>
  static void eliminate_tracing(P& ptc) {
    // FIXME
    // std::size_t id = ptc.template get<pid>();
    // apt::setbits<0, 1>(id, false);
    // ptc.template set<pid>(id);
  }

  template <typename P>
  static bool is_traced(P& ptc) {
    std::size_t id = ptc.template get<pid>();
    return apt::getbits<0, 1>(id);
  }

  template <typename P>
  static bool was_traced(P& ptc) {
    std::size_t id = ptc.template get<pid>();
    id = apt::getbits<1 + Nbits_wr, Nbits_sn>(id);
    return id != 0;
  }

  template <typename R, template <typename> class S>
  static std::unordered_map<species, std::size_t> parse_tracing(
      const map<array<R, S>>& particles);

  template <int DGrid, typename R, template <typename> class S>
  static std::string save_tracing(
      std::string prefix, const int num_parts,
      const std::optional<dye::Ensemble<DGrid>>& ens_opt, int timestep,
      const particle::map<particle::array<R, S>>& particles);

 private:
  std::size_t _wr{};  // world rank
  std::unordered_map<species, std::size_t>
      _counter{};  // NOTE new item guaranteed to be zero-initialized

  static TracingManager& _data() {
    static TracingManager r;
    return r;
  }

  TracingManager() = default;
  TracingManager(const TracingManager&);
  TracingManager(TracingManager&&) noexcept;
  TracingManager& operator=(const TracingManager&);
  TracingManager& operator=(TracingManager&&) noexcept;
  ~TracingManager() = default;
};

}  // namespace particle

namespace particle {
template <int DGrid, typename R, template <typename> class S, typename RJ>
struct Tracer : public Action<DGrid, R, S, RJ> {
 private:
  std::vector<species> _sps;

  bool _is_check_within_range = true;

  pic::Plan _plan{};

  using FMark_t = void (*)(typename array<R, S>::particle_type& ptc,
                           util::Rng<R>& rng);
  FMark_t _marker = nullptr;

 public:
  Tracer* Clone() const override { return new auto(*this); }

  auto& set_marker(FMark_t f) noexcept {
    _marker = f;
    return *this;
  }
  auto& set_species(const std::vector<species>& sps) noexcept {
    _sps = sps;
    return *this;
  }
  auto& set_is_check_within_range(bool a) noexcept {
    _is_check_within_range = a;
    return *this;
  }
  auto& set_plan(const pic::Plan& p) noexcept {
    _plan = p;
    return *this;
  }

  static bool is_within_bounds(
      const typename array<R, S>::particle_type::vec_type& q,
      const apt::array<apt::array<R, 2>, DGrid>& bds) {
    for (int i = 0; i < DGrid; ++i) {
      if (q[i] < bds[i][0] or q[i] >= bds[i][1]) return false;
    }
    return true;
  }

  void operator()(map<array<R, S>>& particles, field::Field<RJ, 3, DGrid>&,
                  std::vector<Particle<R, S>>*, const map<Properties>&,
                  const field::Field<R, 3, DGrid>&,
                  const field::Field<R, 3, DGrid>&,
                  const apt::Grid<R, DGrid>& grid, const dye::Ensemble<DGrid>*,
                  R, int timestep, util::Rng<R>& rng) override {
    if (!_plan.is_do(timestep) or apt::range::is_empty(*this) or !_marker)
      return;

    apt::array<apt::array<R, 2>, DGrid> bds;
    for (int i = 0; i < DGrid; ++i) {
      bds[i][0] = grid[i].absc(apt::range::begin(*this, i));
      bds[i][1] = grid[i].absc(apt::range::end(*this, i));
    }

    for (auto sp : _sps) {
      for (auto ptc : particles[sp]) {  // TODOL semantics
        if (!ptc.is(flag::exist) or
            (_is_check_within_range and !is_within_bounds(ptc.q(), bds)))
          continue;
        _marker(ptc, rng);
      }
    }
  }
};
}  // namespace particle

#endif
