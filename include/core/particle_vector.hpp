#ifndef _PARTICLE_VECTOR_HPP_
#define _PARTICLE_VECTOR_HPP_

#include "core/particle.hpp"
#include <iterator>

namespace particle {
  template < std::size_t Dim_Ptc, typename vector_t, typename T, typename state_t >
  class iterator {
  private:
    vector_t& _vector;
    int _index;

  public:
    using iterator_category = std::random_access_iterator_tag;
    using difference_type = int;
    using value_type = void;
    using reference = Particle< apt::copy_const_t<vector_t, T&>,  Dim_Ptc, apt::copy_const_t<state_t, T&> >;
    using pointer = void;
    iterator( vector_t& vec, int i ) noexcept : _vector(vec), _index(i) {}

    // TODO check _array is the same?
    inline bool operator== ( const iterator& it ) const noexcept {
      return _index == it._index; }

    inline bool operator!= (const iterator& it) const  noexcept {
      return !(*this == it); }

    // prefix ++ // TODO check bounds
    iterator& operator++ () noexcept { ++_index; return *this; }

    iterator& operator+= ( int n ) noexcept { _index += n; return *this; }

    inline reference operator* () const noexcept {
      return reference {
        apt::per_dim::tie<Dim_Ptc> ( [i=_index] ( auto&& x ) { return x[i]; }, _vector._q ),
        apt::per_dim::tie<Dim_Ptc> ( [i=_index] ( auto&& x ) { return x[i]; }, _vector._p ),
        _vector._state[_index] };
    }

  };


  template < typename T, std::size_t Dim_Ptc, typename state_t,
             // ensure only non_cvref qualified types are allowed
             class = std::enable_if_t<
               ( std::is_same_t< T, apt::remove_cvref_t<T> > &&
                 std::is_same_t< state_t, apt::remove_cvref_t<state_t> >
                 ), int> >
  struct vector {
  private:
    std::size_t _capacity = 0;
    std::size_t _size = 0;

    std::array<T*, Dim_Ptc> _q;
    std::array<T*, Dim_Ptc> _p;
    state_t* _state;

  public:
    static constexpr auto DPtc = Dim_Ptc;
    // TODO double this type. Note the & in the end
    using iterator_t = iterator<Dim_Ptc, vector, T&, state_t&>;
    using const_iterator_t = iterator<Dim_Ptc, vector, const T&, const state_t&>;

    vector(std::size_t capacity);
    ~vector();

    inline auto size() const noexcept { return _size; }

    iterator_t begin() noexcept { return iterator_t( *this, 0 ); }
    const_iterator_t begin() const noexcept { return const_iterator_t( *this, 0 ); }

    iterator_t end() noexcept { return iterator_t( *this, size() ); }
    const_iterator_t end() const noexcept { return const_iterator_t( *this, size() ); }

    // TODO check performance
    auto operator[] ( int i ) noexcept { return *( iterator_t( *this, i ) ); }
    auto operator[] ( int i ) const noexcept { return *( const_iterator_t( *this, i ) ); }

    // real particle
    void push_back( const Particle<T, Dim_Ptc, state_t>& ptc );
    void push_back( Particle<T, Dim_Ptc, state_t>&& ptc );

    // virtual particle
    void push_back( const Particle< const T&, Dim_Ptc, const state_t& >& ptc );

  };



}

namespace vec {
  template < typename T, std::size_t DPtc, typename state_t,
             typename PtcRef = vector<T, DPtc, state_t >::interator_t::reference >
  void swap( PtcRef a, PtcRef b ) noexcept;
}


#endif
