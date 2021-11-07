#ifndef _LOGGER_OFSTREAM_HPP_
#define _LOGGER_OFSTREAM_HPP_

#include <fstream>
#include <string>

#include "apt/print.hpp"

namespace lgr {
template <typename CharT = char>
struct ofstream {
 private:
  std::basic_ofstream<CharT> _fout;
  std::string _filename{};
  std::string _indent_old{};
  std::string _indent{};
  bool _is_on = false;

 public:
  inline void set_filename(std::string filename) noexcept {
    _filename = std::move(filename);
  }

  template <typename... Args>
  inline void open(std::string filename, Args&&... args) {
    if (_fout.is_open()) _fout.close();
    _filename = std::move(filename);
    _fout.open(_filename.c_str(), std::forward<Args>(args)...);
    if (_fout.is_open()) _is_on = true;
  }

  template <typename... Args>
  inline void open(Args&&... args) noexcept {
    if (_fout.is_open()) _fout.close();
    _fout.open(_filename.c_str(), std::forward<Args>(args)...);
    if (_fout.is_open()) _is_on = true;
  }

  inline void close() {
    _filename = {};
    if (_fout.is_open()) _fout.close();
    _is_on = false;
  }

  void clear() {
    if (_fout.is_open()) {
      _fout.close();
      _fout.open(_filename.c_str());
    }
  }

  template <typename T>
  ofstream& operator<<(T&& t) {
    if (_is_on) _fout << std::forward<T>(t);
    return *this;
  }

  ofstream& operator<<(std::basic_streambuf<CharT>* sb) {
    if (_is_on) _fout << sb;
    return *this;
  }

  ofstream& operator<<(std::ios_base& (*func)(std::ios_base&)) {
    if (_is_on) _fout << func;
    return *this;
  }

  ofstream& operator<<(std::basic_ios<CharT>& (*func)(std::basic_ios<CharT>&)) {
    if (_is_on) _fout << func;
    return *this;
  }

  ofstream& operator<<(
      std::basic_ostream<CharT>& (*func)(std::basic_ostream<CharT>&)) {
    if (_is_on) _fout << func;
    return *this;
  }

  // start with an indent. % has higher precedence over <<
  // NOTE this assumes t is not std::endl and stuff which will need special
  // treatment as above
  template <typename T>
  ofstream& operator%(T&& t) {
    if (_is_on) _fout << _indent << std::forward<T>(t);
    return *this;
  }

  inline void indent_append(std::string s) {
    _indent_old = _indent;
    _indent += s;
  }

  inline void indent_reset() {
    _indent = _indent_old;
    _indent_old = {};
  }

  inline void turn_on() noexcept { _is_on = true; }
  inline void turn_off() noexcept { _is_on = false; }
};
}  // namespace lgr

#endif
