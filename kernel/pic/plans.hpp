#ifndef _PIC_PLANS_HPP_
#define _PIC_PLANS_HPP_
#include <limits>

namespace pic {
struct Plan {
  bool on = false;
  int start = 0;
  int interval = 100;
  int end = std::numeric_limits<int>::max();

  virtual bool is_do(int timestep) const noexcept {
    return on and timestep >= start and timestep < end and
           ((timestep - start) % interval == 0);
  }
};

Plan sort_ptcs_plan, vitals_plan;

struct ExportPlan : public Plan {
  int num_files = 1;
  int downsample_ratio = 1;
} export_plan;

struct CheckpointPlan : public Plan {
  int num_files = 1;
  int max_num_checkpoints = 1;
  // TODO hourly autosave
  // or ( pic::checkpoint_autosave_hourly &&
  //      autosave.is_save({*pic::checkpoint_autosave_hourly * 3600, "s"}) )
} checkpoint_plan;

struct DynamicLoadBalancePlan : public Plan {
  std::size_t target_load = 100000;
} load_balance_plan;

struct ProfilingPlan : public Plan {
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
} profiling_plan;

struct SaveTracingPlan : public Plan {
  int num_files = 1;
} save_tracing_plan;

int print_timestep_to_stdout_interval = 100;

std::optional<int (*)(int)>
    init_replica_deploy{};  // take in ensemble label and return the intended
                            // number of replicas in that ensemble
}  // namespace pic
#endif
