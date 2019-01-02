#ifndef  _FIELD_HPP_
#define  _FIELD_HPP_

#include <array>
#include <vector>

template < typename T, int DGrid >
struct FieldComponent {
private:
  std::vector<T> data;

public:
  const std::array< bool, DGrid > offset;

  FieldComponent( size_t size, std::array< bool, DGrid > ofs )
    : offset(ofs) {
    data = std::vector<T>( size, T() );
    data.shrink_to_fit();
  }

  inline T& operator[] ( int index ) { return data[index]; }
  inline const T& operator[] ( int index ) const { return data[index]; }

};

template < typename T, int Dim_Field, int Dim_Grid >
struct Field {
  using data_type = T;
  static constexpr int DField = Dim_Field;
  static constexpr int DGrid = Dim_Grid;

  const std::array< int, Dim_Grid > anchor;
  const std::array< int, Dim_Grid > stride;
  std::array < FieldComponent<T, Dim_Grid>, Dim_Field > components;

  inline auto& operator[] ( int i_comp ) { return components[i_comp]; }
  inline const auto& operator[] ( int i_comp ) const { return components[i_comp]; }

};

#endif
