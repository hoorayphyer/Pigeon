#ifndef _PARTICLE_SORTER_HPP_
#define _PARTICLE_SORTER_HPP_
#include <algorithm>

#include "apt/grid.hpp"
#include "particle/array.hpp"

namespace particle {
template <typename T, template <typename> class Spec>
void erase_nonexist(array<T, Spec>& ptcs) {
  int head = 0;
  int tail = ptcs.size() - 1;
  while (head <= tail) {
    if (ptcs[head].is(flag::exist)) {
      ++head;
    } else {
      if (!ptcs[tail].is(flag::exist)) {
        --tail;
      } else {
        ptcs[head] = std::move(ptcs[tail]);
        ++head;
        --tail;
      }
    }
  }
  // when exit head = new size
  ptcs.resize(head);
}

template <typename T, typename F>
void permute(std::vector<T>& perm, F f) {
  // perm has the meaning of the perm[i]-th element before permutation is now at
  // i-th place
  for (T i = 0; i < perm.size(); ++i) {
    // for each i, the following is dealing with a cycle starting at i.
    auto j = i;
    while (perm[j] != i) {
      f(perm[j], j);
      auto k = perm[j];
      perm[j] = j;
      j = k;
    }
    perm[j] = j;
  }
}

template <typename T, int DGrid, template <typename> class Spec>
void sort(array<T, Spec>& ptcs, const apt::Grid<T, DGrid>& grid) {
  erase_nonexist(ptcs);

  std::vector<std::size_t> perm(ptcs.size());
  for (std::size_t i = 0; i < perm.size(); ++i) perm[i] = i;

  constexpr T tile = 1.0;

  auto comp = [&](const std::size_t& a, const std::size_t& b) {
    for (int c = 0; c < Spec<T>::Dim; ++c) {
      auto diff = (ptcs[a].q(c) - ptcs[b].q(c)) / grid[c].delta();
      if (diff < -tile) return true;
      if (diff > tile) return false;
    }
    return false;
  };

  std::sort(perm.begin(), perm.end(), comp);

  auto f = [&ptcs](std::size_t a, std::size_t b) mutable {
    for (int c = 0; c < Spec<T>::Dim; ++c)
      std::swap(ptcs[a].q(c), ptcs[b].q(c));
    for (int c = 0; c < Spec<T>::Dim; ++c)
      std::swap(ptcs[a].p(c), ptcs[b].p(c));
    std::swap(ptcs[a].frac(), ptcs[b].frac());
    std::swap(ptcs[a].state(), ptcs[b].state());
  };

  permute(perm, f);
}
}  // namespace particle

#endif
