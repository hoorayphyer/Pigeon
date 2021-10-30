#include "silopp/silo_optlist.hpp"

#include <silo.h>

#include <unordered_map>

namespace silo {
enum class Val_t : char { FLOAT = 0, DOUBLE, INT, INT3 };

struct OptMetadata {
  int raw_opt;
  Val_t type;
};

std::unordered_map<Opt, OptMetadata> meta = []() {
  std::unordered_map<Opt, OptMetadata> res;
  res[Opt::TIME] = {DBOPT_TIME, Val_t::FLOAT};
  res[Opt::DTIME] = {DBOPT_DTIME, Val_t::DOUBLE};
  res[Opt::CYCLE] = {DBOPT_CYCLE, Val_t::INT};
  res[Opt::LO_OFFSET] = {DBOPT_LO_OFFSET, Val_t::INT3};
  res[Opt::HI_OFFSET] = {DBOPT_HI_OFFSET, Val_t::INT3};
  res[Opt::BASEINDEX] = {DBOPT_BASEINDEX, Val_t::INT3};

  return res;
}();
}  // namespace silo

namespace silo {
OptVal::OptVal(Opt id) {
  switch (meta.at(id).type) {
    case Val_t::FLOAT:
      _numeric = (float)0;
      break;
    case Val_t::DOUBLE:
      _numeric = (double)0;
      break;
    case Val_t::INT3:
      _int3.reset(new int[3]);
      for (int i = 0; i < 3; ++i) _int3[i] = 0;
      break;
    default:
      _numeric = (int)0;
  }
}
};  // namespace silo

namespace silo {
void optlist_free(DBoptlist* optlist) { DBFreeOptlist(optlist); }

OptList::OptList() { _p.reset(DBMakeOptlist(3)); }

OptList::OptList(const OptList& other)
    : std::unordered_map<Opt, OptVal>(
          static_cast<const std::unordered_map<Opt, OptVal>&>(other)) {
  _p.reset(DBMakeOptlist(3));
}

OptList::operator DBoptlist*() const {
  for (auto& [id, val] : *this) {
    const auto& [raw_opt, val_type] = meta.at(id);
    if (val_type == Val_t::INT3) {
      DBAddOption(_p.get(), raw_opt, val._int3.get());
    } else {
      std::visit([opt = raw_opt, ptr = _p.get()](
                     auto& v) { DBAddOption(ptr, opt, (void*)&v); },
                 val._numeric);
    }
  }

  return _p.get();
}
}  // namespace silo
