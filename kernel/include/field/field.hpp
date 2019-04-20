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

    inline T& operator() ( const apt::Index<Mesh_t::NDim>& i_bulk ) {
      return _data [_mesh.linearized_index_of_whole_mesh(i_bulk) ];
    }

    inline T operator() ( const apt::Index<Mesh_t::NDim>& i_bulk ) const {
      return _data [_mesh.linearized_index_of_whole_mesh(i_bulk) ];
    }

    inline T& operator[] ( int i ) {
      return _data[i];
    }

    inline T operator[] ( int i ) const {
      return _data[i];
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
    Mesh<DGrid> _mesh;
    apt::array< apt::array< offset_t, DGrid >, DField > _offset;

  public:
    using element_type = T;
    static constexpr int NDim = DField;

    Field() = default;

    Field( const Mesh<DGrid>& mesh ) : _mesh(mesh) {
      int size = 1;
      for ( int i = 0; i < DGrid; ++i ) size *= _mesh.extent()[i];
      apt::foreach<0,DField>
        ( [size] ( auto& comp ) {
            comp.reserve( size );
            comp.resize( size );
            for ( auto& elm : comp ) elm = static_cast<T>(0);
          }, _comps );
    }

    inline Field& set_offset( int component, const apt::array< offset_t, DGrid >& offset ) noexcept {
      _offset[component] = offset;
      return *this;
    }

    inline const auto& mesh() const noexcept { return _mesh;}

    inline const auto operator[] ( int i ) const noexcept {
      return Component< const std::vector<T>, Mesh<DGrid>, apt::array<offset_t, DGrid> >
        ( _comps[i], _mesh, _offset[i] );
    }

    inline auto operator[] ( int i ) noexcept {
      return Component< std::vector<T>, Mesh<DGrid>, apt::array<offset_t, DGrid> >
        ( _comps[i], _mesh, _offset[i] );
    }

    inline void reset() {
      apt::foreach<0, DField>
        ( []( auto& comp ) {
            for ( auto& elm : comp ) elm = 0.0;
          }, _comps );
    }


  };
}

#endif
