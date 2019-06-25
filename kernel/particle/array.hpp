#ifndef _PARTICLE_ARRAY_HPP_
#define _PARTICLE_ARRAY_HPP_

#include "particle/virtual_particle.hpp"
#include "particle/particle.hpp"
#include <iterator>
#include <vector>

namespace ckpt {
  template < typename T, template < typename > class PtcSpecs >
  struct ParticleArrayCkpt;
}

namespace particle {
  template < typename T, template < typename > class Specs >
  struct array {
  private:
    apt::array<std::vector<T>, Specs<T>::Dim> _q;
    apt::array<std::vector<T>, Specs<T>::Dim> _p;
    std::vector<typename Specs<T>::state_type> _state;

  public:
    using value_type = T;
    static constexpr int DPtc = Specs<T>::Dim;

    using particle_type = vParticle< T, Specs >;

    using const_particle_type = vParticle< const T, Specs >;

    friend class ckpt::ParticleArrayCkpt<T,Specs>;

    template < bool isConst >
    class iterator {
    private:
      using array_t = std::conditional_t<isConst, const array, array>;
      array_t& _array;
      int _index;

      // ref_tuple is not default constructible. So recursion is used instead
      template < typename Vec, int I = 0 >
      constexpr auto build_ref_tuple( Vec& x, int i, apt::ref_tuple<decltype(x[0][0]),I> res = {} ) noexcept {
        if constexpr ( I == DPtc ) return res;
        else return build_ref_tuple<Vec,I+1>( x, i, std::tuple_cat( std::move(res), std::tie(x[I][i]) ) );
      }

    public:
      using iterator_category = std::random_access_iterator_tag;
      using difference_type = int;
      using value_type = void;
      using reference = std::conditional_t<isConst, array::const_particle_type, array::particle_type >;
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
        return reference( build_ref_tuple(_array._q, _index),
                          build_ref_tuple(_array._p, _index),
                          _array._state[_index] );
      }

    };

    friend class iterator<false>;
    friend class iterator<true>;

    inline auto size() const noexcept { return _state.size(); }

    auto begin() noexcept { return iterator<false>( *this, 0 ); }
    auto begin() const noexcept { return iterator<true>( *this, 0 ); }

    auto end() noexcept { return iterator<false>( *this, size() ); }
    auto end() const noexcept { return iterator<true>( *this, size() ); }

    // TODO check performance
    auto operator[] ( int i ) noexcept {
      return *( iterator<false>( *this, i ) );
    }
    auto operator[] ( int i ) const noexcept {
      return *( iterator<true>( *this, i ) );
    }

    template < typename Ptc >
    void push_back( Ptc&& ptc ) {
      auto f = []( auto& arr, auto&& x ) {
                 arr.push_back(std::forward<decltype(x)>(x)); };
      apt::foreach<0,DPtc> ( f, _q, ptc.q() );
      apt::foreach<0,DPtc> ( f, _p, ptc.p() );
      f( _state, ptc.state() );
    }

    template < typename... Args >
    inline void emplace_back( Args&&... args ) {
      push_back( Particle<T,Specs>(std::forward<Args>(args)...) );
    }

    // TODO emplace_back with braced-lists

    // NOTE from is inclusive, to is exclusive. from can be larger than to.
    void erase( unsigned int from, unsigned int to );

    void resize(std::size_t size);
  };

}

namespace std {
  template < typename T, template <typename> class Specs  >
  class back_insert_iterator<particle::array<T,Specs>> {
  private:
    particle::array<T,Specs>& _arr;
    int _index;

    using self_type = back_insert_iterator<particle::array<T,Specs>>;
  public:
    using iterator_category = std::output_iterator_tag;
    using difference_type = void;
    using value_type = void;
    using reference = void;
    using pointer = void;

    explicit back_insert_iterator( particle::array<T,Specs>& arr ) noexcept : _arr(arr) {}

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
