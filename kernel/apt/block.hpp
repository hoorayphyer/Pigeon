#ifndef  _APT_BLOCK_HPP_
#define  _APT_BLOCK_HPP_

#include "apt/index.hpp"
#include <iterator>

namespace apt {
  // A wrapper over int so as to utilize type-safe checks
  struct BlockIteratorEnd {
  private:
    const int _end = 0;

  public:
    constexpr BlockIteratorEnd( int end ) noexcept : _end(end) {}
    constexpr operator int() const noexcept { return _end; }
  };

  template < int D >
  struct BlockIterator {
  private:
    Index<D> _ijk{};
    Index<D> _extent;

  public:
    using difference_type = void;
    using value_type = void;
    using reference = const Index<D>&;
    using pointer = void;
    using iterator_category = std::forward_iterator_tag;

    constexpr BlockIterator( const Index<D>& extent ) noexcept
      : _extent( extent ) {}

    constexpr bool operator!= ( const BlockIteratorEnd& end ) const noexcept {
      return _ijk[D-1] != static_cast<int>(end);
    }

    // NOTE separating ++ ijk and ijk %= _extent is the key to make iteration stoppable
    constexpr BlockIterator& operator++() noexcept {
      ++ (_ijk[0]);
      if constexpr ( D > 1 ) {
          _ijk[0] %= _extent[0];
          if ( 0 != _ijk[0] ) return *this;
          else ++ (_ijk[1]);
        }
      if constexpr ( D > 2 ) {
          _ijk[1] %= _extent[1];
          if ( 0 != _ijk[1]) return *this;
          else ++ (_ijk[2]);
        }
      return *this;
    }

    constexpr BlockIterator operator++(int) noexcept {
      auto res = *this;
      ++(*this);
      return res;
    }

    constexpr reference operator*() noexcept { return _ijk; }
  };
}

// NOTE to be used in range-based for, Block and BlockIterator has to be separated classes because there is auto __begin = __range.begin(), which prevents __range.begin() from returning itself as iterator
namespace apt {
  template < int D >
  struct Block {
  private:
    apt::Index<D> _extent{};

  public:
    constexpr Block( const apt::Index<D>& extent ) noexcept : _extent(extent) {}

    constexpr auto begin() const noexcept { return BlockIterator<D>(_extent);}

    constexpr BlockIteratorEnd end() const noexcept {
      // deal with empty or invalid block
      for ( int i = 0; i < D; ++i ) {
        if ( _extent[i] < 1 ) return {0};
      }
      return {_extent[D-1]};
    }

  };
}


#endif
