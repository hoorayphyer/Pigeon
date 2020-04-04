#ifndef _PIC_MODULES_
#define _PIC_MODULES_
#include <limits>

namespace pic {
  struct ModuleBase {
    bool on = false;
    int start = 0;
    int end = std::numeric_limits<int>::max();
    int interval = 100;

    virtual bool is_do( int timestep ) const noexcept {
      return on and timestep >= start and timestep < end and (timestep % interval == 0 );
    }
  };

  ModuleBase mod_sort_ptcs, mod_vitals;

  struct ExportModule : public ModuleBase {
    int num_files = 1;
    int downsample_ratio = 1;
  } mod_export;

  struct CheckpointModule : public ModuleBase {
    int num_files = 1;
    int max_num_checkpoints = 1;
    // TODO hourly autosave
    // or ( pic::checkpoint_autosave_hourly &&
    //      autosave.is_save({*pic::checkpoint_autosave_hourly * 3600, "s"}) )
  } mod_checkpoint;

  struct DynamicLoadBalanceModule : public ModuleBase {
    std::size_t target_load = 100000;
  } mod_load_balance;

  struct ProfilingModule : public ModuleBase {
    std::optional<int> max_entries {100};
    bool (*is_qualified) () = [](){return true;};

    bool is_reached_max_entries( int timestep ) const noexcept {
      if ( !max_entries ) return false;
      static int last = timestep;
      if ( timestep >= last + *max_entries * interval ) {
        last = timestep;
        return true;
      }
      return false;
    }
  } mod_profiling;

  struct TracingModule : public ModuleBase {
    int num_files = 1;
  } mod_tracing;

  int print_timestep_to_stdout_interval = 100;

  std::optional<int (*) ( int )> init_replica_deploy {}; // take in ensemble label and return the intended number of replicas in that ensemble
}
#endif
