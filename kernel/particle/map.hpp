#pragma once

#include <array>
#include <optional>
#include <stdexcept>
#include <vector>

#include "particle/species_predef.hpp"

namespace particle {
// NOTE this data structure map is to achieve the following
// - the order of keys is always fixed by the enum values
// - values are put in contiguous memory so as to ease MPI communication
template <typename Val>
struct map {
 private:
  std::vector<Val> _data;     // ordered by enum values
  std::vector<species> _sps;  // same as above
  std::array<std::optional<int>, NUM_SPECIES>
      _ind;  // index is enum value cast to T
  using T = std::underlying_type_t<species>;

 public:
  inline bool has(species sp) const noexcept {
    return static_cast<bool>(_ind[static_cast<T>(sp)]);
  }

  template <typename... Args>
  void insert(species sp, Args&&... args) {
    if (has(sp)) {
      _data[*_ind[static_cast<T>(sp)]] = {std::forward<Args>(args)...};
    } else {
      _data.emplace_back(std::forward<Args>(args)...);
      _sps.emplace_back(sp);
      int idx = 0;
      for (int i = 0; i < static_cast<T>(sp); ++i) {
        if (_ind[i]) ++idx;
      }

      _ind[static_cast<T>(sp)] = idx;

      for (int i = static_cast<T>(sp) + 1; i < NUM_SPECIES; ++i) {
        if (_ind[i]) {
          // first swap _data elements. then update _ind[i]
          std::swap(_data[*_ind[i]], _data.back());
          std::swap(_sps[*_ind[i]], _sps.back());
          ++(*_ind[i]);
        }
      }
    }
  }

  // this overload is to enable brace-initialization
  inline void insert(species sp, const Val& v) { insert<const Val&>(sp, v); }

  void erase(species sp) {
    if (!has(sp)) return;
    int idx = *_ind[static_cast<T>(sp)];
    // swap with later species one by one
    for (int i = idx + 1; i < NUM_SPECIES; ++i) {
      if (_ind[i]) {
        // first swap _data elements. then update _ind[i]
        std::swap(_data[*_ind[i]], _data[idx]);
        std::swap(_sps[*_ind[i]], _sps[idx]);
        std::swap(*_ind[i], idx);
      }
    }
    _data.pop_back();
    _sps.pop_back();
    _ind[static_cast<T>(sp)].reset();
  }

  inline Val& operator[](species sp) {
    if (!has(sp)) {
      throw std::runtime_error("Species " +
                               std::to_string(static_cast<int>(sp)) +
                               " doesn't exist in this map!");
    }
    return _data[*_ind[static_cast<T>(sp)]];
  }

  inline const Val& operator[](species sp) const {
    if (!has(sp)) {
      throw std::runtime_error("Species " +
                               std::to_string(static_cast<int>(sp)) +
                               " doesn't exist in this map!");
    }
    return _data[*_ind[static_cast<T>(sp)]];
  }

  inline auto& data() { return _data; }
  inline const auto& data() const { return _data; }

  inline const auto begin() const noexcept { return _sps.cbegin(); }

  inline const auto end() const noexcept { return _sps.cend(); }

  inline auto size() const noexcept { return _data.size(); }
};
}  // namespace particle
