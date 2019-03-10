#ifndef  _FIELD_HPP_
#define  _FIELD_HPP_

#include "field/offset.hpp"
#include "field/mesh.hpp"
#include <array>
#include <vector>

namespace field {
  template < typename Data_t, typename Mesh_t, typename Offset_t, typename T = decltype(Data_t()[0])  >
  struct Component {
  private:
    Data_t& _data;
    const Mesh_t& _mesh;
    const Offset_t& _offset;

  public:
    constexpr Component( Data_t& data, const Mesh_t& mesh, const Offset_t& offset ) noexcept
      : _data(data), _mesh(mesh), _offset(offset) {}

    template < typename... Indices >
    inline T& operator() ( const Indices&... indices_bulk  ) {
      return _data [_mesh.linearized_index_of_whole_mesh(indices_bulk...) ];
    }

    template < typename... Indices >
    inline T operator() ( const Indices&... indices_bulk  ) const {
      return _data [_mesh.linearized_index_of_whole_mesh(indices_bulk...) ];
    }

    inline auto offset() const noexcept { return _offset; }
  };

}

namespace field {

  template < typename T, int DField, int DGrid >
  struct Field {
  private:
    std::array < std::vector<T>, DField > _comps;
    Mesh<T,DGrid> _mesh;
    std::array< std::array< offset_t, DGrid >, DField > _offset;

  public:
    using element_type = T;
    static constexpr int NDim = DField;

    Field( const Mesh<T,NDim>& mesh ) : _mesh(mesh) {
      for( auto& comp : _comps ) {
        comp.resize( _mesh.size() );
        for ( auto& elm : comp ) elm = static_cast<T>(0);
        comp.shrink_to_fit();
      }
    }

    template < int Comp >
    Field& set_offset( const std::array< offset_t, DGrid >& offset ) {
      std::get<Comp>(_offset) = offset;
      return *this;
    }

    inline const auto& mesh() const { return _mesh;}

    template < int Comp >
    inline const auto c() const noexcept {
      return Component< const std::vector<T>, Mesh<T,NDim>, std::array<offset_t, DGrid> >
        ( std::get<Comp>(_comps), _mesh, std::get<Comp>(_offset) );
    }

    template < int Comp >
    inline auto c() noexcept {
      return Component< std::vector<T>, Mesh<T,NDim>, std::array<offset_t, DGrid> >
        ( std::get<Comp>(_comps), _mesh, std::get<Comp>(_offset) );
    }

  };
}

namespace std {
  // define this so as to be used in apt::foreach
  template < int I, typename T, int DF, int DG >
  constexpr decltype(auto) get ( const field::Field<T,DF,DG>& f ) noexcept {
    return f.template c<I>();
  }

  template < int I, typename T, int DField, int DGrid >
  constexpr decltype(auto) get ( field::Field<T,DField,DGrid>& f ) noexcept {
    return f.template c<I>();
  }
}

#endif
