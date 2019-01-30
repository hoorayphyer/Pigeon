#ifndef  _FIELD_HPP_
#define  _FIELD_HPP_

#include <array>
#include <vector>

template < typename T, int DGrid >
struct FieldComponent {
private:
  std::vector<T> _data;

public:
  const std::array< T, DGrid > offset;

  FieldComponent( size_t size, std::array< bool, DGrid > ofs )
    : offset(ofs) {
    _data = std::vector<T>( size, T() );
    _data.shrink_to_fit();
  }

  inline T& operator[] ( int index ) { return _data[index]; }
  inline const T& operator[] ( int index ) const { return _data[index]; }

};

template < typename T, int Dim_Field, int Dim_Grid >
struct Field {
  using data_type = T;
  static constexpr auto DField = Dim_Field;
  static constexpr auto DGrid = Dim_Grid;

  const std::array< int, Dim_Grid > anchor;
  const std::array< int, Dim_Grid > stride;
  std::array < FieldComponent<T, Dim_Grid>, Dim_Field > components;

  // inline auto& operator[] ( int i_comp ) { return components[i_comp]; }
  // inline const auto& operator[] ( int i_comp ) const { return components[i_comp]; }

  //------ int... inds_local -------
  // TODO indices
  template < int Comp >
  inline decltype(auto) c( int... inds ) {
    static_assert( sizeof...(inds) >= Field::DGrid, "Not enough indices provided" );
    // TODO inner_product
    auto&& ind_linear = tum::inner_product( field.stride, std::forward_as_tuple(std::move(inds)...) );
    return std::get<Comp>(components)[std::move(ind_linear)];
  }

  template < int Comp >
  inline decltype(auto) c( int... inds ) const {
    static_assert( sizeof...(inds) >= Field::DGrid, "Not enough indices provided" );
    // TODO inner_product
    auto&& ind_linear = tum::inner_product( field.stride, std::forward_as_tuple(std::move(inds)...) );
    return std::get<Comp>(components)[std::move(ind_linear)];
  }

  //------ int... inds_global -------
  template < int Comp >
  inline decltype(auto) c( std::array<int,DGrid> I_global ) {
    // TODO inner_product
    auto&& ind_linear = tum::inner_product( field.stride, std::move(I_global) - anchor );
    return std::get<Comp>(components)[std::move(ind_linear)];
  }

  template < int Comp >
  inline decltype(auto) c( std::array<int,DGrid> I_global ) const {
    // TODO inner_product
    auto&& ind_linear = tum::inner_product( field.stride, std::move(I_global) - anchor );
    return std::get<Comp>(components)[std::move(ind_linear)];
  }


};

template < typename T, int DField, int DGrid >
Field< T, DField, DGrid >
make_field( std::array< int, DGrid > anchor,
            std::array< int, DGrid > extent,
            std::array< std::array< bool, DGrid >, DField > offsets );


#endif
