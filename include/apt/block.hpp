#ifndef  _APT_BLOCK_HPP_
#define  _APT_BLOCK_HPP_

#include "apt/index.hpp"
#include <iterator>

namespace apt {
  template < int D >
  struct BlockIterator {
  private:
    Index<D> _ijk{};
    const Index<D>& _extent;

  public:
    using difference_type = void;
    using value_type = void;
    using reference = const Index<D>&;
    using pointer = void;
    using iterator_category = std::forward_iterator_tag;

    constexpr BlockIterator( const Index<D>& extent ) noexcept
      : _extent( extent ) {}

    constexpr bool operator!= ( const Index<D>& idx ) const noexcept {
      return _ijk != idx;
    }

    constexpr BlockIterator& operator++() noexcept {
      ( ++ (_ijk[0]) ) %= _extent[0];
      if constexpr ( D > 1 ) {
          if ( 0 != _ijk[0] ) return *this;
          else ( ++ (_ijk[1]) ) %= _extent[1];
        }
      if constexpr ( D > 2 ) {
          if ( 0 != _ijk[1]) return *this;
          else ( ++ (_ijk[2]) ) %= _extent[2];
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

namespace apt {
  template < int D >
  struct Block {
  private:
    apt::Index<D> _extent{};

  public:
    constexpr Block( apt::Index<D> extent ) noexcept : _extent(std::move(extent)) {}

    constexpr auto begin() const noexcept { return BlockIterator<D>(_extent);}

    constexpr auto end() const noexcept {
      auto res = _extent;
      apt::foreach<1,D> ( [](auto& x){ --x; }, res );
      return res;
    }

  };
}


#endif
