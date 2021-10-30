#ifndef _DYE_SCATTER_LOAD_
#define _DYE_SCATTER_LOAD_

#include "apt/pair.hpp"

namespace dye {
constexpr apt::pair<int> scatter_load(int total_load, int this_receiver_idx,
                                      int num_receivers) noexcept {
  int from = 0, num = 0;

  int share = total_load / num_receivers;
  int remain = total_load % num_receivers;

  int adj_idx = (this_receiver_idx + num_receivers) % num_receivers;

  if (adj_idx < remain) {
    num = share + 1;
    from = adj_idx * (share + 1);
  } else {
    num = share;
    from = adj_idx * share + remain;
  }
  return {from, num};
}

}  // namespace dye
#endif
