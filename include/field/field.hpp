#ifndef  _FIELD_HPP_
#define  _FIELD_HPP_

#include "field/offset.hpp"
#include "field/mesh.hpp"
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

    inline const auto& data() const noexcept { return _data; }
    inline auto& data() noexcept { return _data; }
  };

}

namespace field {

  template < typename T, int DField, int DGrid >
  struct Field {
  private:
    apt::array < std::vector<T>, DField > _comps;
    Mesh<T,DGrid> _mesh;
    apt::array< apt::array< offset_t, DGrid >, DField > _offset;

  public:
    using element_type = T;
    static constexpr int NDim = DField;

    Field( const Mesh<T,DGrid>& mesh ) : _mesh(mesh) {
      for( auto& comp : _comps ) {
        comp.resize( _mesh.size() );
        for ( auto& elm : comp ) elm = static_cast<T>(0);
        comp.shrink_to_fit();
      }
    }

    inline Field& set_offset( int component, const apt::array< offset_t, DGrid >& offset ) noexcept {
      _offset[component] = offset;
    }

    inline const auto& mesh() const { return _mesh;}

    inline const auto operator[] ( int i ) const noexcept {
      return Component< const std::vector<T>, Mesh<T,DGrid>, apt::array<offset_t, DGrid> >
        ( _comps[i], _mesh, _offset[i] );
    }

    inline auto operator[] ( int i ) noexcept {
      return Component< std::vector<T>, Mesh<T,DGrid>, apt::array<offset_t, DGrid> >
        ( _comps[i], _mesh, _offset[i] );
    }

  };
}

#endif
