#ifndef _PARTICLE_STATE_HPP_
#define _PARTICLE_STATE_HPP_

#include "particle/species_predef.hpp"
#include "particle/flags_predef.hpp"
#include "apt/ternary.hpp"
#include "apt/integer_class.hpp"
#include <type_traits>
#include <tuple>

namespace particle {
  using destination = apt::Integer<int,0>;
  using birthplace = apt::Integer<unsigned int,0>;
  using serial_number= apt::Integer<unsigned int,1>;
}

namespace particle {
  // convention: Pos counts from the right. (Pos+N, Pos].
  // NOTE bit shift operator << returns of type of the left operand. So if LHS has fewer bits than RHS, the return result will simply be zero, i.e., information is lost, hence all the casts
  template < std::size_t Pos, std::size_t N, typename U, typename T >
  constexpr U getbits( const T& x ) noexcept {
    return static_cast<U>( ( x >> Pos ) & ~( ~static_cast<std::size_t>(0) << N ) );
  }

  template < std::size_t Pos, std::size_t N,
             typename U, typename T >
  constexpr void setbits( T& x, U y ) noexcept {
    // NOTE assume y has all nonzero bits in the last N bits
    x &= ~( ~( ~static_cast<std::size_t>(0) << N ) << Pos );
    x |= ( ( static_cast<std::size_t>(y) & ~( ~static_cast<std::size_t>(0) << N ) ) << Pos );
  }

  struct layout {
  private:
    using ordering = std::tuple<species, destination, flag, birthplace, serial_number>;
    static constexpr auto sizing = std::make_tuple(3, 5, 16, 16, 24); // TODO remove the hard coded sum to 64

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
  private:
    template < typename Attr, typename Arg >
    inline auto getbits( Arg&& x) const noexcept {
      return particle::getbits< layout::begin<Attr>(), layout::size<Attr>(), std::underlying_type_t<Attr> >( std::forward<Arg>(x) );
    }

    template < typename Attr, typename... Args >
    inline void setbits( Args&&... x) noexcept {
      particle::setbits< layout::begin<Attr>(), layout::size<Attr>() >( std::forward<Args>(x)... );
    }


    inline void set_impl( const flag& attr ) noexcept {
      using U = std::underlying_type_t<flag>;
      state() |= ( static_cast<std::size_t>(1) << (static_cast<U>(attr) + layout::begin<flag>() ) );
    }

    // The goal is to have an integer-to-ternary-array converter such that int(0) means {C,C,C,....}.
    inline void set_impl( destination attr ) noexcept {
      // set the sign bit, 0 for >= 0, 1 for < 0
      particle::setbits< layout::begin<destination>() + layout::size<destination>() - 1, 1 >( state(), bool( attr < 0 ) );
      attr = ( attr < 0 ) ? -attr : attr;
      // set the rest
      particle::setbits< layout::begin<destination>(), layout::size<destination>() - 1 >( state(), attr );
    }

    template < typename Attr >
    inline void set_impl( const Attr& attr ) noexcept {
      setbits< Attr >( state(), static_cast<std::underlying_type_t<Attr>>(attr) );
    }

    template < typename Attr >
    inline Attr get_impl() const noexcept {
      if constexpr ( std::is_same_v<Attr, destination> ) {
        using U = std::underlying_type_t<destination>;
        U x = particle::getbits< layout::begin<destination>(), layout::size<destination>() - 1, U >( state() );
        x *= ( particle::getbits< layout::begin<destination>() + layout::size<destination>() - 1, 1, U >( state() ) ) ? -1 : 1;
       return {x};
}
      else
        return static_cast<Attr>( getbits<Attr>(state()) );
    }

  public:
    using self_t = StateExpression<Ptc, T>;
    template < int >
    friend class migrInt;

    constexpr T& state() noexcept { return static_cast<Ptc&>(*this).state(); }
    constexpr T state() const noexcept { return static_cast<const Ptc&>(*this).state(); }

    template < typename Attr, typename... Attrs >
    inline std::enable_if_t< !std::is_same_v<Attr, destination>, self_t&> // disable user interface to set destination
    set( const Attr& attr, const Attrs&... attrs ) noexcept {
      set_impl(attr);
      if constexpr( sizeof...(Attrs) > 0 ) set(attrs...);
      return *this;
    }

    inline self_t& reset( const flag& fl ) noexcept {
      using U = std::underlying_type_t<flag>;
      // NOTE must cast 1 to std::size_t to avoid being wiped out
      state() &= ~( static_cast<std::size_t>(1) << (static_cast<U>(fl) + layout::begin<flag>() ) );
      return *this;
    }

    template < typename Attr >
    inline self_t& reset() noexcept {
      set_impl<Attr>(Attr{});
      return *this;
    }

    template < typename Attr >
    inline std::enable_if_t<!std::is_same_v<Attr, destination>, Attr>
    get() const noexcept {
      return get_impl<Attr>();
    }

    template < typename Attr >
    inline bool is( const Attr& attr ) const noexcept {
      if constexpr ( std::is_same_v<Attr, flag> ) {
        using U = std::underlying_type_t<flag>;
        return ( state() >> layout::begin<flag>() ) & ( static_cast<U>(1) << static_cast<U>(attr) );
      }
      else
        return get<Attr>() == attr;
    }


  };

}

#endif
