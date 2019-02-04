#ifndef _PARTICLE_ARRAY_HPP_
#define _PARTICLE_ARRAY_HPP_

#include "particle/particle.hpp"
#include <iterator>

namespace particle {
  template < std::size_t Dim_Ptc, typename array_t, typename T >
  class iterator {
  private:
    array_t& _array;
    int _index;

  public:
    using iterator_category = std::random_access_iterator_tag;
    using difference_type = int;
    using value_type = void;
    using reference = Particle< apt::copy_const_t<array_t, T&>,  Dim_Ptc >;
    using pointer = void;
    iterator( array_t& arr, int i ) noexcept : _array(arr), _index(i) {}

    // TODO check _array is the same?
    inline bool operator== ( const iterator& it ) const noexcept {
      return _index == it._index; }

    inline bool operator!= (const iterator& it) const  noexcept {
      return !(*this == it); }

    // prefix ++ // TODO check bounds
    iterator& operator++ () noexcept { ++_index; return *this; }

    iterator& operator+= ( int n ) noexcept { _index += n; return *this; }

    inline reference operator* () const noexcept {
      // TODO expression template
      // return reference {
        // apt::per_dim::tie<Dim_Ptc> ( [i=_index] ( auto&& x ) { return x[i]; }, _array._q ),
        // apt::per_dim::tie<Dim_Ptc> ( [i=_index] ( auto&& x ) { return x[i]; }, _array._p ),
        // _array._state[_index] };
    }

  };


  template < typename T, std::size_t Dim_Ptc >
  struct array {
    static_assert( std::is_same_v< T, apt::remove_cvref_t<T> >, "only non_cvref qualified types are allowed" );
  private:
    std::size_t _capacity = 0;
    std::size_t _size = 0;

    std::array<T*, Dim_Ptc> _q;
    std::array<T*, Dim_Ptc> _p;
    state_t<T>* _state;

  public:
    static constexpr auto DPtc = Dim_Ptc;
    // TODO double this type. Note the & in the end
    using iterator_t = iterator<Dim_Ptc, array, T& >;
    using const_iterator_t = iterator<Dim_Ptc, array, const T&>;

    array(std::size_t capacity);
    ~array();

    inline auto size() const noexcept { return _size; }

    iterator_t begin() noexcept { return iterator_t( *this, 0 ); }
    const_iterator_t begin() const noexcept { return const_iterator_t( *this, 0 ); }

    iterator_t end() noexcept { return iterator_t( *this, size() ); }
    const_iterator_t end() const noexcept { return const_iterator_t( *this, size() ); }

    // TODO check performance
    auto operator[] ( int i ) noexcept { return *( iterator_t( *this, i ) ); }
    auto operator[] ( int i ) const noexcept { return *( const_iterator_t( *this, i ) ); }

    // real particle
    void push_back( const Particle<T, DPtc>& ptc ); // NOTE cannot use Dim_Ptc here!!!
    void push_back( Particle<T, DPtc>&& ptc );

    // virtual particle TODO move virtual particle??
    void push_back( const Particle< const T&, DPtc >& ptc );

  };



}

namespace apt {
  template < typename T, std::size_t DPtc >
  void swap( typename particle::array<T, DPtc >::interator_t::reference a,
             typename particle::array<T, DPtc >::interator_t::reference b ) noexcept;
}

namespace std {
  template < typename T, std::size_t DPtc >
  class back_insert_iterator<particle::array<T,DPtc>> {
  private:
    particle::array<T,DPtc>& _c;
    int _index;

  public:
    using iterator_category = std::output_iterator_tag;
    using difference_type = void;
    using value_type = void;
    using reference = void;
    using pointer = void;

    explicit back_insert_iterator( particle::array<T,DPtc>& c ) noexcept : _c(c) {}

    auto& operator= ( const particle::Particle<T, DPtc>& ptc ) { _c.push_back(ptc); return *this; }
    auto& operator= ( const particle::Particle< const T&, DPtc >& ptc ) { _c.push_back(ptc); return *this; }
    auto& operator= ( particle::Particle<T, DPtc >&& ptc ) { _c.push_back( std::move(ptc) ); return *this; }

    inline auto& operator++ () noexcept { return *this; }
    inline auto& operator++ (int) noexcept { return *this; }
    inline auto& operator* () noexcept { return *this; }

  };

}


#endif
