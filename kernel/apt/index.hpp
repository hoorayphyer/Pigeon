#ifndef  _APT_INDEX_HPP_
#define  _APT_INDEX_HPP_

#include "apt/array.hpp"
#include "apt/foreach.hpp"

namespace apt {
  template < int D >
  using Index = array<int,D>;
}

template < int D >
constexpr apt::Index<D> operator+( const apt::Index<D>& ind_a, const apt::Index<D>& ind_b ) noexcept {
  apt::Index<D> res;
  apt::foreach<0,D>
    ( []( auto& r, auto a, auto b ) { r = a + b; }, res, ind_a, ind_b );
  return res;
}

template < int D >
constexpr apt::Index<D> operator-( const apt::Index<D>& ind_a, const apt::Index<D>& ind_b ) noexcept {
  apt::Index<D> res;
  apt::foreach<0,D>
    ( []( auto& r, auto a, auto b ) { r = a - b; }, res, ind_a, ind_b );
  return res;
}

template < int D >
constexpr apt::Index<D>& operator+=( apt::Index<D>& ind_a, const apt::Index<D>& ind_b ) noexcept {
  apt::foreach<0,D>
    ( []( auto& a, auto b ) { a += b; }, ind_a, ind_b );
  return ind_a;
}

template < int D >
constexpr bool operator==( const apt::Index<D>& ind_a, const apt::Index<D>& ind_b ) noexcept {
  bool res = true;
  apt::foreach<0,D>
    ( [&res]( auto a, auto b ) { res = ( res && (a == b) ); }, ind_a, ind_b );
  return res;
}

template < int D >
constexpr bool operator!=( const apt::Index<D>& ind_a, const apt::Index<D>& ind_b ) noexcept {
  return !( ind_a == ind_b );
}

namespace apt {
  struct Longidx { // longitudianl index
  private:
    const int _dir = 0;
    int _val = 0;

  public:
    Longidx() = default;
    Longidx(const Longidx&) = default;
    Longidx(Longidx&&) noexcept = default;

    constexpr Longidx& operator=( const Longidx& ) = default;
    constexpr Longidx& operator=( Longidx&& ) noexcept = default;

    constexpr Longidx( int dir, int val = 0) noexcept : _dir(dir), _val(val) {}
    constexpr Longidx& operator=( int v ) noexcept { _val = v; return *this; }
    constexpr operator int() const noexcept { return _val; }

    constexpr const int dir() const noexcept { return _dir; }

    template < int D >
    constexpr Index<D> operator+( Index<D> x ) const noexcept {
      x[_dir] += _val; return x;
    }

    constexpr Longidx& operator++() noexcept { ++_val; return *this; }
    constexpr Longidx& operator--() noexcept { --_val; return *this; }
    constexpr Longidx operator++(int) noexcept { auto res = *this; ++_val; return res; }
    constexpr Longidx operator--(int) noexcept { auto res = *this; --_val; return res; }

    constexpr bool operator<(int a) const noexcept { return _val < a; }
    constexpr bool operator<=(int a) const noexcept { return _val <= a; }
    constexpr bool operator>(int a) const noexcept { return _val > a; }
    constexpr bool operator>=(int a) const noexcept { return _val >= a; }
    constexpr bool operator!=(int a) const noexcept { return _val != a; }
    constexpr bool operator==(int a) const noexcept { return _val == a; }
  };
}

template < int D >
constexpr apt::Index<D> operator+( const apt::Index<D>& ind, const apt::Longidx& l ) noexcept { return l + ind;}

template < int D >
constexpr apt::Index<D> operator-( apt::Index<D> ind, const apt::Longidx& l ) noexcept {
  ind[l.dir()] -= l;
  return ind;
}



#endif
