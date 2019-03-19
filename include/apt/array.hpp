#ifndef  _APT_ARRAY_HPP_
#define  _APT_ARRAY_HPP_

namespace apt {
  template < typename T, int D >
  struct array {
    T _data[D] {};

    using element_type = T;
    static constexpr int NDim = D;

    // TODOL bound checks.
    constexpr T operator[] ( int i ) const noexcept { return _data[i]; }
    constexpr T& operator[] ( int i ) noexcept { return _data[i]; }

    bool operator== ( const array& other ) const noexcept {
      bool res = true;
      for ( int i = 0; i < D; ++i )
        res &= ( _data[i] == other[i] );
      return res;
    }

    constexpr T* data() noexcept { return _data;}
    constexpr const T* data() const noexcept { return _data;}

    constexpr int size() const noexcept { return NDim; }
  };

  // a C-array is not allowed to have size 0
  template < typename T >
  struct array<T,0> {
    static constexpr int NDim = 0;
  };
}

#endif
