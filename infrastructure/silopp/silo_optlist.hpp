#ifndef _SILO_OPTLIST_HPP_
#define _SILO_OPTLIST_HPP_

#include <memory>
#include <unordered_map>
#include <variant>
#include <vector>

typedef struct DBoptlist_ DBoptlist;
using RawList = DBoptlist;

namespace silo {
enum class Opt : int {
  TIME = 0,
  DTIME,
  CYCLE,
  LO_OFFSET,
  HI_OFFSET,
  BASEINDEX
};
}

namespace silo {
void optlist_free(RawList*);

struct OptList;

struct OptVal {
 private:
  std::variant<int, float, double> _numeric;
  std::unique_ptr<int[]>
      _int3;  // NOTE variant is not allowed to hold reference, array, or void

 public:
  friend struct OptList;
  OptVal(Opt opt_id);

  OptVal() = default;
  OptVal(const OptVal& other) {
    _numeric = other._numeric;
    if (other._int3) {
      if (!_int3) _int3.reset(new int[3]);
      for (int i = 0; i < 3; ++i) _int3[i] = other._int3[i];
    }
  }
  OptVal(OptVal&&) = default;
  ~OptVal() = default;

  // TODO does this change the type of val? Seems no from reference, which IS
  // what we want
  template <typename T>
  std::enable_if_t<!std::is_same_v<T, std::vector<int>>, OptVal&> operator=(
      const T& t) noexcept {
    std::visit([&t](auto& v) { v = t; }, _numeric);
    return *this;
  }

  OptVal& operator=(const std::vector<int>& t) {
    for (int i = 0; i < t.size(); ++i) _int3[i] = t[i];
    return *this;
  }
};

struct OptList : public std::unordered_map<Opt, OptVal> {
 private:
  // used only to auto manage resources when converting to RawList*. No need to
  // be copied in copy constructor
  std::unique_ptr<RawList, void (*)(RawList*)> _p{nullptr, optlist_free};

 public:
  operator RawList*() const;

  OptVal& operator[](Opt opt_id) {
    emplace(opt_id, OptVal(opt_id));
    return at(opt_id);
  }

  OptList();
  OptList(const OptList& other);
  OptList(OptList&&) = default;
  ~OptList() = default;
};

}  // namespace silo

#endif
