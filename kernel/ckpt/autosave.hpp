#ifndef _CHECKPOINT_AUTOSAVE_HPP_
#define _CHECKPOINT_AUTOSAVE_HPP_

#include <optional>
#include <queue>
#include <string>

#include "timer/timer.hpp"

namespace ckpt {
struct Autosave {
 private:
  tmr::Timestamp _stamp;
  std::queue<std::string> _q;

 public:
  bool is_save(tmr::Duration interval) {
    bool res = false;
    if (_stamp.lapse() > interval) res = true;
    _stamp = {};
    return res;
  }

  std::optional<std::string> add_checkpoint(std::string new_ckpt_dir,
                                            std::size_t max_ckpts) {
    std::optional<std::string> res;
    if (_q.size() >= max_ckpts) {
      res.emplace(_q.front());
      _q.pop();
    }
    _q.push(new_ckpt_dir);
    return res;
  }
};
}  // namespace ckpt

#endif
