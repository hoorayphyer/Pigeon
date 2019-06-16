#ifndef _MODULE_RANGE_
#define _MODULE_RANGE_
namespace pic {
  struct ModuleRange {
    bool is_on {};
    int init_ts = 0;
    int interval = 100;
  };
}
#endif
