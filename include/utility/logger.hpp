#ifndef  _LOGGER_HPP_
#define  _LOGGER_HPP_

#include <iostream>
#include <fstream>
#include <string>





namespace Logger {
  extern std::ofstream debug;
  extern bool isDebugActiveRank; // true if one wants the rank to output debug
  extern std::string debugLogFile;
  extern std::string dest; // either "file" or "screen"
  template <typename T>
  void debug_print(const T& t) {
    if ( !isDebugActiveRank ) return;
    if ("file" == dest)
      debug << t << std::endl;
    else if ( "screen" == dest )
      std::cout << t << std::endl;
    else
      std::cerr << "unknown debug destination " << dest << std::endl;
  }

  template <typename First, typename... Rest>
  void debug_print(const First& first, const Rest&... rest) {
    if ( !isDebugActiveRank ) return;
    if ("file" == dest)
      debug << first << " ";
    else if ( "screen" == dest )
      std::cout << first << " ";
    else
      std::cerr << "unknown debug destination " << dest << std::endl;
    debug_print(rest...); // recursive call using pack expansion syntax
  }

  void debug_clear_log(); // delete the debug log so far

}

extern std::string debug_indent;
// for doing debug_indent_str * 2
std::string operator*( std::string s, int n );


#endif
