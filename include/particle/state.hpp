#ifndef _PARTICLE_STATE_HPP_
#define _PARTICLE_STATE_HPP_

#include "particle/traits.hpp"
#include <type_traits>
#include <tuple>

namespace particle {
  namespace trace {
    struct birthplace {
    private:
      unsigned int _value;
    public:
      birthplace( unsigned int a) noexcept : _value(a) {}
      operator unsigned int() noexcept { return _value; }
    };

    struct serial_number {
    private:
      unsigned int _value;
    public:
      serial_number( unsigned int a) noexcept : _value(a) {}
      operator unsigned int() noexcept { return _value; }
    };
  }
}


namespace particle {

  template < typename T >
  struct state {
    static_assert( 8 * sizeof( std::remove_reference_t<T> ) >= 64 );
  private:
    struct layout {
    private:
      using ordering = std::tuple<species, flag, trace::birthplace, trace::serial_number>;
      static constexpr auto sizing = std::make_tuple(2, 16, 16, 30);

    public:
      template < typename Attr , std::size_t I = 0 >
      static constexpr std::size_t begin() noexcept {
        static_assert ( (I < std::tuple_size_v<ordering>), "unknown Attr" );
        using elm_t = std::tuple_element_t<I, ordering>;
        if constexpr( std::is_same_v< Attr, elm_t > ) return 0;
        else return size<elm_t> + begin<Attr, I+1>();
      }

      template < typename Attr , std::size_t I = 0 >
      static constexpr std::size_t size() noexcept {
        static_assert ( (I < std::tuple_size_v<ordering>), "unknown Attr" );
        using elm_t = std::tuple_element_t<I, ordering>;
        if constexpr( std::is_same_v< Attr, elm_t > ) return std::get<I>(sizing);
        else return size<Attr, I+1>();
      }
    };

    // convention: Pos counts from the right. (Pos+N, Pos].
    template < std::size_t Pos, std::size_t N, typename T,
               class = std::enable_if_t< std::is_unsigned_v<T>, int > >
    constexpr auto getbits( const T& x ) noexcept {
      return ( x >> Pos ) & ~( ~0 << N );
    }

    template < std::size_t Pos, std::size_t N,
               typename T, typename U,
               class = std::enable_if_t<
                 std::is_unsigned_v<T> &&
                 std::is_unsigned_v<U> , int > >
    constexpr void setbits( T& x, U y ) noexcept {
      x &= ~( ~( ~0 << N ) << Pos );
      x |= ( ( y &= ~( ~0 << N ) ) << Pos );
    }

  protected:
    T _state; // NOTE use T here so when T is nonref, _state will hold the data

  public:
    state() noexcept = default;

    template < typename Attr, typename... Attrs >
    state( const Attr& attr, const Attrs&... attrs ) noexcept {
      set(attr);
      if constexpr( sizeof...(Attrs) > 0 ) return set(attrs...);
    }


    template < typename Attr >
    inline auto& set( const Attr& a ) noexcept {
      setbits<layout::begin<Attr>(), layout::size<Attr>() >( _state, static_cast<std::underlying_type_t<Attr>(a)> );
      return *this;
    }

    template <>
    inline auto& set<flag>( const flag& fl ) noexcept {
      using utt = std::underlying_type_t<flag>;
      _state |= ( 1 << (static_cast<utt>(fl) + layout::begin<flag>() ) );
      return *this;
    }

    inline auto& reset( const flag& fl ) noexcept {
      using utt = std::underlying_type_t<flag>;
      _state &= ~( 1 << (static_cast<utt>(fl) + layout::begin<flag>() ) );
      return *this;
    }

    template < typename Attr >
    inline auto get() const noexcept {
      return static_cast<Attr>( getbits< layout::begin<Attr>(), layout::size<Attr>() >(_state) );
    }

    template < typename Attr >
    inline bool is( const Attr& a ) const noexcept {
      return get<Attr>() == a;
    }

    template <>
    inline bool is<flag>( const flag& fl ) const noexcept {
      return _state & ( 1 << (static_cast<utt>(fl) + layout::begin<flag>()) );
    }

  };

}

#endif
