#ifndef  _APT_BLOCK_HPP_
#define  _APT_BLOCK_HPP_

#include "apt/index.hpp"
#include <iterator>

namespace apt {
  template <int> struct Block;

  // Specialiation of longitudianl block
  template <>
  struct Block<1> {
  private:
    int _begin = 0;
    int _end = 0;
    Longidx _i {};

  public:
    constexpr Block( int longi, int begin, int end ) noexcept
      : _begin(std::move(begin)), _end(std::move(end)), _i(longi,_begin) {}

    constexpr Block begin() const noexcept { return {_i.dir(),_begin,_end}; }

    constexpr int end() const noexcept {
      if ( _end <= _begin ) return _begin;
      else return _end;
    }

    using difference_type = void;
    using value_type = void;
    using reference = Longidx;
    using pointer = void;
    using iterator_category = std::forward_iterator_tag;

    constexpr bool operator!= ( int end ) const noexcept { return _i != end; }

    constexpr Block& operator++() noexcept { ++_i; return *this; }

    constexpr Block operator++(int) noexcept {auto res = *this; ++(*this); return res;}

    constexpr reference operator*() noexcept { return _i; }
  };
}

namespace apt {
  template < int D >
  struct Block {
  private:
    Index<D> _begin{};
    Index<D> _end{};
    Index<D> _ijk{};

  public:
    using end_type = int;
    constexpr Block( apt::Index<D> begin, apt::Index<D> end ) noexcept
      : _begin(std::move(begin)), _end(std::move(end)), _ijk(_begin) {}

    constexpr Block begin() const noexcept { return {_begin,_end};}

    constexpr end_type end() const noexcept {
      // deal with empty or invalid block
      for ( int i = 0; i < D; ++i ) {
        if ( _end[i] <= _begin[i] ) return {_begin[D-1]};
      }
      return {_end[D-1]};
    }

    using difference_type = void;
    using value_type = void;
    using reference = const Index<D>&;
    using pointer = void;
    using iterator_category = std::forward_iterator_tag;

    constexpr bool operator!= ( const end_type& end ) const noexcept {
      return _ijk[D-1] != end;
    }

    // NOTE separating ++ ijk and ijk %= _extent is the key to make iteration stoppable
    constexpr Block& operator++() noexcept {
      ++ (_ijk[0]);
      if constexpr ( D > 1 ) {
          if ( _ijk[0] != _end[0] ) return *this;
          else {
            _ijk[0] = _begin[0];
            ++ (_ijk[1]);
          }
        }
      if constexpr ( D > 2 ) {
          if ( _ijk[1] != _end[1] ) return *this;
          else {
            _ijk[1] = _begin[1];
            ++ (_ijk[2]);
          }
        }
      return *this;
    }

    constexpr Block operator++(int) noexcept {
      auto res = *this;
      ++(*this);
      return res;
    }

    constexpr reference operator*() noexcept { return _ijk; }
  };

  template< int D >
  constexpr apt::Block<D> project_out(int longi, apt::Index<D> b, apt::Index<D> e ) noexcept {
    b[longi] = 0;
    e[longi] = 1;
    return {b,e};
  }
}

#endif
