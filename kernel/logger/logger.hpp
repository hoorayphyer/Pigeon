#ifndef  _LOGGER_HPP_
#define  _LOGGER_HPP_

#include "apt/print.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

namespace lgr {
  template < typename CharT = char >
  struct ofstream {
  private:
    std::basic_ofstream<CharT> _fout;
    std::string indent_old{};
    std::string indent{};
    bool _is_on = true;

  public:
    template < typename... Args >
    inline auto open( Args&&... args ) {
      return _fout.open( std::forward<Args>(args)... );
    }

    inline auto close() { return _fout.close(); }

    template < typename T >
    ofstream& operator<< ( T&& t ) {
      if ( _is_on ) _fout << std::forward<T>(t);
      return *this;
    }

    ofstream& operator<<( std::basic_streambuf<CharT>* sb) {
      if ( _is_on ) _fout << sb;
      return *this;
    }

    ofstream& operator<<( std::ios_base& (*func)(std::ios_base&) ) {
      if ( _is_on ) _fout << func;
      return *this;
    }

    ofstream& operator<<( std::basic_ios<CharT>& (*func)(std::basic_ios<CharT>&) ) {
      if ( _is_on ) _fout << func;
      return *this;
    }

    ofstream& operator<<( std::basic_ostream<CharT>& (*func)( std::basic_ostream<CharT>& ) ) {
      if ( _is_on ) _fout << func;
      return *this;
    }

    // start with an indent. % has higher precedence over <<
    // NOTE this assumes t is not std::endl and stuff which will need special treatment as above
    template < typename T >
    ofstream& operator% ( T&& t ) {
      if ( _is_on ) _fout << indent << std::forward<T>(t);
      return *this;
    }

    inline void indent_append( std::string s ) {
      indent_old = indent;
      indent += s;
    }

    inline void indent_reset() {
      indent = indent_old;
      indent_old = {};
    }

    inline void turn_on() noexcept { _is_on = true; }
    inline void turn_off() noexcept { _is_on = false; }

  };
}

namespace lgr {
  extern ofstream<> file;
  extern ofstream<> debug;

  extern std::ostream& out;
  extern std::ostream& err;
  extern std::ostringstream str;
}

#endif
