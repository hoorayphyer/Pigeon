#ifndef _MODULE_RANGE_
#define _MODULE_RANGE_
namespace pic {
  struct ModuleRange {
    bool is_on {};
    int init_ts = 0;
    int interval = 100;

    inline bool is_do( int timestep ) const noexcept {
      return is_on && timestep >= init_ts && (timestep % interval == 0 );
    }
  };

}
#endif
