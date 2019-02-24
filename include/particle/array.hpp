#ifndef _PARTICLE_ARRAY_HPP_
#define _PARTICLE_ARRAY_HPP_

#include "particle/virtual_particle.hpp"
#include <iterator>
#include <vector>

namespace particle {
  template < typename array_t >
  class iterator {
  private:
    array_t& _array;
    int _index;

  public:
    using iterator_category = std::random_access_iterator_tag;
    using difference_type = int;
    using value_type = void;
    using reference = vParticle< typename array_t::value_type, array_t::DPtc, typename array_t::state_type >;
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

    inline reference operator* () noexcept {
      // TODO generalize to all DPtc
      using vVec = apt::vVec<typename array_t::value_type, array_t::DPtc>;
      auto f = [i=_index] ( auto& x ) {
                 return vVec( std::get<0>(x)[i], std::get<1>(x)[i], std::get<2>(x)[i]);
               };
      return reference( f( _array._q ),
                        f( _array._p ),
                        _array._state[_index] );
    }

  };

  template < typename T, int Dim_Ptc, typename state_t >
  struct array {
  private:
    std::array<std::vector<T>, Dim_Ptc> _q;
    std::array<std::vector<T>, Dim_Ptc> _p;
    std::vector<state_t> _state;

  public:
    using value_type = T;
    static constexpr auto DPtc = Dim_Ptc;
    using state_type = state_t;

    friend class iterator< array >;
    friend class iterator< const array >;

    inline auto size() const noexcept { return _state.size(); }

    auto begin() noexcept { return iterator( *this, 0 ); }
    auto begin() const noexcept { return iterator( *this, 0 ); }

    auto end() noexcept { return iterator( *this, size() ); }
    auto end() const noexcept { return iterator( *this, size() ); }

    // TODO check performance
    auto operator[] ( int i ) noexcept {
      return *( iterator( *this, i ) );
    }
    auto operator[] ( int i ) const noexcept {
      return *( iterator( *this, i ) );
    }

    template < typename Ptc >
    void push_back( const PtcExpression<Ptc>& ptc );
    template < typename Ptc >
    void push_back( PtcExpression<Ptc>&& ptc );

    // NOTE from is inclusive, to is exclusive. from can be larger than to.
    void erase( int from, int to );

    void resize(std::size_t size);

  };



}

namespace std {
  template < typename T, int DPtc, typename state_t >
  class back_insert_iterator<particle::array<T,DPtc, state_t>> {
  private:
    particle::array<T,DPtc,state_t>& _arr;
    int _index;

    using self_type = back_insert_iterator<particle::array<T,DPtc, state_t>>;
  public:
    using iterator_category = std::output_iterator_tag;
    using difference_type = void;
    using value_type = void;
    using reference = void;
    using pointer = void;

    explicit back_insert_iterator( particle::array<T,DPtc, state_t>& arr ) noexcept : _arr(arr) {}

    template < typename Ptc >
    self_type& operator= ( Ptc&& ptc ) {
      _arr.push_back(std::forward<Ptc>(ptc));
      return *this;
    }

    inline self_type& operator++ () noexcept { return *this; }
    inline self_type& operator++ (int) noexcept { return *this; }
    inline self_type& operator* () noexcept { return *this; }

  };

}


#endif
