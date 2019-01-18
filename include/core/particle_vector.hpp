#ifndef _PARTICLE_VECTOR_HPP_
#define _PARTICLE_VECTOR_HPP_

#include "core/particle.hpp"
#include <vector>
#include <iterator>

namespace particle {
  template < std::size_t Dim_Ptc, typename vector_t, typename... Ts  >
  class iterator {
  private:
    vector_t& _vector;
    int _index;

  public:
    using iterator_category = std::random_access_iterator_tag;
    using difference_type = int;
    // using value_type = Particle< Dim_Ptc, Ts... >;
    using value_type = void;
    using reference = Particle< Dim_Ptc, std::conditional_t<
                                           std::is_const_v<vector_t>,
                                           const Ts&..., Ts&...> >;
    using pointer = void;

    iterator( vector_t& vec, int i ) noexcept : _vector(vec), _index(i) {}

    // TODO check _array is the same?
    inline bool operator== ( const iterator& it ) const noexcept {
      return _index == it._index; }

    inline bool operator!= (const iterator& it) const  noexcept {
      return !(*this == it); }

    // prefix ++ // TODO check bounds
    iterator& operator++ () { ++_index; return *this; }

    iterator& operator+= ( int n ) { _index += n; return *this; }

    inline reference operator* () const {
      return reference {
        vec::per_dim::tie<Dim_Ptc> ( [i=_index] ( auto&& x ) { return x[i]; }, _vector.q ),
        vec::per_dim::tie<Dim_Ptc> ( [i=_index] ( auto&& x ) { return x[i]; }, _vector.p ),
        _vector.state[_index] };
    }

  };

  template < std::size_t Dim_Ptc, typename T >
  struct vector {
  private:
    using encoded_bits_t = unsigned long long; // need at least 64 bits
    std::array<std::vector<T>, Dim_Ptc> q;
    std::array<std::vector<T>, Dim_Ptc> p;
    std::vector<encoded_bits_t> state;

  public:
    static constexpr auto DPtc = Dim_Ptc;
    // TODO double this type. Note the & in the end
    using iterator_t = iterator<Dim_Ptc, vector, T&, encoded_bits_t&>;
    using const_iterator_t = iterator<Dim_Ptc, vector, const T&, const encoded_bits_t&>;

    inline auto size() const noexcept { return state.size(); }

    iterator_t begin() { return iterator_t( *this, 0 ); }
    const_iterator_t begin() const { return const_iterator_t( *this, 0 ); }

    iterator_t end() { return iterator_t( *this, size() ); }
    const_iterator_t end() const { return const_iterator_t( *this, size() ); }

    auto operator[] ( int i ) {
      return *( iterator_t( *this, i ) );
    }

    auto operator[] ( int i ) const {
      return *( const_iterator_t( *this, i ) );
    }

  };
}

#endif
