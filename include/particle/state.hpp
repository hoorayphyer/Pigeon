#ifndef _PARTICLE_STATE_HPP_
#define _PARTICLE_STATE_HPP_

#include "particle/species_predef.hpp"
#include "particle/flags_predef.hpp"
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
  // convention: Pos counts from the right. (Pos+N, Pos].
  template < std::size_t Pos, std::size_t N, typename T >
  constexpr auto getbits( const T& x ) noexcept {
    static_assert(std::is_unsigned_v<T>);
    return ( x >> Pos ) & ~( ~0 << N );
  }

  template < std::size_t Pos, std::size_t N,
             typename T, typename U >
  constexpr void setbits( T& x, U y ) noexcept {
    static_assert(std::is_unsigned_v<T>);
    static_assert(std::is_unsigned_v<U>);
    x &= ~( ~( ~0 << N ) << Pos );
    x |= ( ( y &= ~( ~0 << N ) ) << Pos );
  }

  struct layout {
  private:
    using ordering = std::tuple<species, flag, trace::birthplace, trace::serial_number>;
    static constexpr auto sizing = std::make_tuple(2, 16, 16, 30);

  public:
    template < typename Attr, std::size_t I = 0 >
    static constexpr std::size_t size() noexcept {
      static_assert ( (I < std::tuple_size_v<ordering>), "unknown Attr" );
      using elm_t = std::tuple_element_t<I, ordering>;
      if constexpr( std::is_same_v< Attr, elm_t > ) return std::get<I>(sizing);
      else return size<Attr, I+1>();
    }

    template < typename Attr, std::size_t I = 0 >
    static constexpr std::size_t begin() noexcept {
      static_assert ( (I < std::tuple_size_v<ordering>), "unknown Attr" );
      using elm_t = std::tuple_element_t<I, ordering>;
      if constexpr( std::is_same_v< Attr, elm_t > ) return 0;
      else return size<elm_t>() + begin<Attr, I+1>();
    }
  };

  template < typename Ptc, typename T >
  struct StateExpression {
    using self_t = StateExpression<Ptc, T>;

    constexpr T& state() noexcept { return static_cast<Ptc&>(*this).state(); }
    constexpr T state() const noexcept { return static_cast<const Ptc&>(*this).state(); }

    template < typename Attr, typename... Attrs >
    inline self_t& set( const Attr& attr, const Attrs&... attrs ) noexcept {

      if constexpr ( std::is_same_v<Attr, flag> ) {
      state() |= ( 1 << (static_cast<std::underlying_type_t<Attr>>(attr) + layout::begin<flag>() ) );
      } else if ( std::is_same_v<Attr, species>) {
        setbits< layout::begin<Attr>(), layout::size<Attr>() >( state(), static_cast<std::underlying_type_t<Attr>>(attr) );
      } else {
        setbits< layout::begin<Attr>(), layout::size<Attr>() >( state(), (unsigned int)attr );
      }

      if constexpr( sizeof...(Attrs) > 0 ) set(attrs...);
      return *this;
    }

    inline self_t& reset( const flag& fl ) noexcept {
      state() &= ~( 1 << (static_cast<std::underlying_type_t<flag>>(fl) + layout::begin<flag>() ) );
      return *this;
    }

    template < typename Attr >
    inline Attr get() const noexcept {
      return static_cast<Attr>( getbits< layout::begin<Attr>(), layout::size<Attr>() >(state()) );
    }

    template < typename Attr >
    inline bool is( const Attr& attr ) const noexcept {
      if constexpr ( std::is_same_v<Attr, flag> )
        return state() & ( 1 << (static_cast<std::underlying_type_t<flag>>(attr) + layout::begin<flag>()) );
      else
        return get<Attr>() == attr;
    }


  };


}

#endif
