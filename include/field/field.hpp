#ifndef  _FIELD_HPP_
#define  _FIELD_HPP_

#include <array>
#include <vector>

namespace field {
  struct offset_t {
    bool val;
    constexpr operator bool() noexcept { return val; }
    template < typename T >
    constexpr operator T() noexcept { return static_cast<T>( 0.5 * val ); }
  };

  constexpr offset_t INSITU{false};
  constexpr offset_t MIDWAY{true};

  template < typename T, int DGrid >
  struct Component : std::vector<T> {
    const std::array< offset_t, DGrid > offset;
  };

}

namespace field {
  template < typename T, int Dim_Field, int Dim_Grid >
  struct Field : public std::array < Component<T, Dim_Grid>, Dim_Field > {
    using data_type = T;
    static constexpr auto DField = Dim_Field;
    static constexpr auto DGrid = Dim_Grid;

    const std::array< int, Dim_Grid > anchor;
    const std::array< int, Dim_Grid > extent;
    const std::array< int, Dim_Grid + 1 > stride;

    //------ int... inds_global -------
    template < int Comp >
    inline auto& c( const std::array<int,DGrid>& I_global ) {
      auto&& ind_linear = std::get<0>(stride) * ( std::get<0>(I_global) - std::get<0>(anchor) );
      if constexpr ( DGrid > 1 )
                     ind_linear += std::get<1>(stride) * ( std::get<1>(I_global) - std::get<1>(anchor) );
      if constexpr ( DGrid > 2 )
                     ind_linear += std::get<2>(stride) * ( std::get<2>(I_global) - std::get<2>(anchor) );
      static_assert( DGrid < 4 );
      return std::get<Comp>(*this)[std::move(ind_linear)];
    }

    template < int Comp >
    inline auto c( const std::array<int,DGrid>& I_global ) const {
      auto&& ind_linear = std::get<0>(stride) * ( std::get<0>(I_global) - std::get<0>(anchor) );
      if constexpr ( DGrid > 1 )
                     ind_linear += std::get<1>(stride) * ( std::get<1>(I_global) - std::get<1>(anchor) );
      if constexpr ( DGrid > 2 )
                     ind_linear += std::get<2>(stride) * ( std::get<2>(I_global) - std::get<2>(anchor) );
      static_assert( DGrid < 4 );
      return std::get<Comp>(*this)[std::move(ind_linear)];
    }

    //------ local index
    template < int Comp >
    inline auto& c( int i_local_linear ) {
      return std::get<Comp>(*this)[i_local_linear];
    }

    template < int Comp >
    inline auto c( int i_local_linear ) const {
      return std::get<Comp>(*this)[i_local_linear];
    }

  };
}

#endif
