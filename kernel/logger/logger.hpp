#ifndef  _LOGGER_HPP_
#define  _LOGGER_HPP_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

namespace lgr {
  struct ofstream : public std::ofstream {
  private:
    std::string indent_old{};
    std::string indent{};
  public:
    // start with an indent. % has higher precedence over <<
    template < typename T >
    std::ofstream& operator% ( T&& t ) {
      *this << indent << std::forward<T>(t);
      return static_cast<std::ofstream&>(*this);
    }

    inline void indent_append( std::string s ) {
      indent_old = indent;
      indent += s;
    }

    inline void indent_reset() {
      indent = indent_old;
      indent_old = {};
    }

  };
}

namespace lgr {
  extern ofstream file;
  extern ofstream debug;

  auto& out = std::cout;
  auto& err = std::cerr;
  extern std::ostringstream str;
}

#endif
