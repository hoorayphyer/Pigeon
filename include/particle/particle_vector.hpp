#ifndef _PARTICLE_VECTOR_HPP_
#define _PARTICLE_VECTOR_HPP_

#include "particle/particle.hpp"
#include <iterator>

namespace particle {
  template < std::size_t Dim_Ptc, typename vector_t, typename T >
  class iterator {
  private:
    vector_t& _vector;
    int _index;

  public:
    using iterator_category = std::random_access_iterator_tag;
    using difference_type = int;
    using value_type = void;
    using reference = Particle< apt::copy_const_t<vector_t, T&>,  Dim_Ptc >;
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
      // TODO expression template
      // return reference {
        // apt::per_dim::tie<Dim_Ptc> ( [i=_index] ( auto&& x ) { return x[i]; }, _vector._q ),
        // apt::per_dim::tie<Dim_Ptc> ( [i=_index] ( auto&& x ) { return x[i]; }, _vector._p ),
        // _vector._state[_index] };
    }

  };


  template < typename T, std::size_t Dim_Ptc,
             // ensure only non_cvref qualified types are allowed
             class = std::enable_if_t< std::is_same_v< T, apt::remove_cvref_t<T> >, int> >
  struct vector {
  private:
    std::size_t _capacity = 0;
    std::size_t _size = 0;

    std::array<T*, Dim_Ptc> _q;
    std::array<T*, Dim_Ptc> _p;
    state_t<T>* _state;

  public:
    static constexpr auto DPtc = Dim_Ptc;
    // TODO double this type. Note the & in the end
    using iterator_t = iterator<Dim_Ptc, vector, T& >;
    using const_iterator_t = iterator<Dim_Ptc, vector, const T&>;

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
    void push_back( const Particle<T, Dim_Ptc>& ptc );
    void push_back( Particle<T, Dim_Ptc>&& ptc );

    // virtual particle
    void push_back( const Particle< const T&, Dim_Ptc >& ptc );

  };



}

namespace apt {
  template < typename T, std::size_t DPtc >
  void swap( typename particle::vector<T, DPtc >::interator_t::reference a,
             typename particle::vector<T, DPtc >::interator_t::reference b ) noexcept;
}

namespace std {
  template < typename T, std::size_t DPtc >
  class back_insert_iterator<particle::vector<T,DPtc>> {
  private:
    particle::vector<T,DPtc>& _c;
    int _index;

  public:
    using iterator_category = std::output_iterator_tag;
    using difference_type = void;
    using value_type = void;
    using reference = void;
    using pointer = void;

    explicit back_insert_iterator( particle::vector<T,DPtc>& c ) noexcept : _c(c) {}

    template < typename U,
               class = std::enable_if_t< std::is_same_v< T, apt::remove_cvref_t<U> >, int> >
    back_insert_iterator& operator= ( const Particle<U, DPtc >& ptc );

    template < typename U,
               class = std::enable_if_t< std::is_same_v< T, apt::remove_cvref_t<U> >, int> >
    back_insert_iterator& operator= ( Particle<U, DPtc >&& ptc );

    inline auto& operator++ () noexcept { return *this; }
    inline auto& operator++ (int) noexcept { return *this; }
    inline auto& operator* () noexcept { return *this; }

  };

}


#endif
