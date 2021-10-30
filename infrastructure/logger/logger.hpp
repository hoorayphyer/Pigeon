#ifndef _LOGGER_HPP_
#define _LOGGER_HPP_

#include <iostream>
#include <sstream>

#include "logger/ofstream.hpp"

namespace lgr {
extern ofstream<> file;
extern ofstream<> debug;

extern std::ostream& out;
extern std::ostream& err;
extern std::ostringstream str;
}  // namespace lgr

#endif
