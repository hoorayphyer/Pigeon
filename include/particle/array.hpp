#ifndef _PARTICLE_ARRAY_HPP_
#define _PARTICLE_ARRAY_HPP_

#include "particle/virtual_particle.hpp"
#include "apt/vec_from_function.hpp"
#include <iterator>
#include <memory>

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
      auto f =
        [i=_index] ( auto& x ) -> typename array_t::value_type& // force returning by lref
        { return x[i]; };
      return reference( apt::make_vff<array_t::DPtc>( f, _array._q ),
                        apt::make_vff<array_t::DPtc>( f, _array._p ),
                        _array._state[_index] );
    }

  };

  template < typename T, int Dim_Ptc, typename state_t >
  struct array {
  private:
    std::size_t _capacity = 0;
    std::size_t _size = 0;

    std::array<std::unique_ptr<T[]>, Dim_Ptc> _q;
    std::array<std::unique_ptr<T[]>, Dim_Ptc> _p;
    std::unique_ptr<state_t[]> _state;

  public:
    using value_type = T;
    static constexpr auto DPtc = Dim_Ptc;
    using state_type = state_t;

    friend class iterator< array >;
    friend class iterator< const array >;

    array(std::size_t capacity);
    ~array();

    inline auto size() const noexcept { return _size; }

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

    // TODO change them to ParticleExpression
    template < typename Ptc >
    void push_back( const PtcExpression<Ptc>& ptc );

    template < typename Ptc >
    void push_back( PtcExpression<Ptc>&& ptc );

    // NOTE from is inclusive, to is exclusive. from can be larger than to.
    void erase( int from, int to );

    // TODO
    void resize(std::size_t size);

  };



}

namespace std {
  template < typename T, int DPtc, typename state_t >
  class back_insert_iterator<particle::array<T,DPtc, state_t>> {
  private:
    particle::array<T,DPtc,state_t>& _c;
    int _index;

  public:
    using iterator_category = std::output_iterator_tag;
    using difference_type = void;
    using value_type = void;
    using reference = void;
    using pointer = void;

    explicit back_insert_iterator( particle::array<T,DPtc, state_t>& c ) noexcept : _c(c) {}

    template < typename Ptc >
    auto& operator= ( Ptc&& ptc ) {
      // use emplace to accommodate for PtcExpression
      _c.emplace_back(std::forward<Ptc>(ptc));
      return *this;
    }

    inline auto& operator++ () noexcept { return *this; }
    inline auto& operator++ (int) noexcept { return *this; }
    inline auto& operator* () noexcept { return *this; }

  };

}


#endif
