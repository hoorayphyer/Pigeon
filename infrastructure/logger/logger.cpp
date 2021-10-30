#include "logger/logger.hpp"

namespace lgr {
std::string indent{};

ofstream<> file;
ofstream<> debug;

std::ostream& out = std::cout;
std::ostream& err = std::cerr;

std::ostringstream str;
}  // namespace lgr
