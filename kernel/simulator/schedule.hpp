#pragma once

#include <limits>
#include <optional>

namespace pic {
struct Schedule {
  bool on = false;
  int start = 0;
  int interval = 100;
  int end = std::numeric_limits<int>::max();

  virtual bool is_do(int timestep) const noexcept {
    return on and timestep >= start and timestep < end and
           ((timestep - start) % interval == 0);
  }
};

struct ExportSchedule : public Schedule {
  int num_files = 1;
};

struct CheckpointSchedule : public Schedule {
  int num_files = 1;
  int max_num_checkpoints = 1;
  // TODO hourly autosave
  // or ( pic::checkpoint_autosave_hourly &&
  //      autosave.is_save({*pic::checkpoint_autosave_hourly * 3600, "s"}) )
};

struct DynamicLoadBalanceSchedule : public Schedule {
  std::size_t target_load = 100000;
};

struct ProfilingSchedule : public Schedule {
  std::optional<int> max_entries{100};
  bool (*is_qualified)() = []() { return true; };

  bool is_reached_max_entries(int timestep) const noexcept {
    if (!max_entries) return false;
    static int last = timestep;
    if (timestep >= last + *max_entries * interval) {
      last = timestep;
      return true;
    }
    return false;
  }
};
}  // namespace pic
