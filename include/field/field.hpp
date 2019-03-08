#ifndef  _FIELD_HPP_
#define  _FIELD_HPP_

#include "field/offset.hpp"
#include "field/mesh.hpp"
#include <array>
#include <vector>

namespace field {
  template < typename T, int DGrid >
  struct Component : std::vector<T> {
    std::array< offset_t, DGrid > offset;
  };

}

namespace field {
  template < typename T, int Dim_Field, int Dim_Grid >
  struct Field {
  private:
    std::array < Component<T, Dim_Grid>, Dim_Field > _comps;
    Mesh<T,Dim_Grid> _mesh;

  public:
    using element_type = T;
    static constexpr int NDim = Dim_Field;

    Field( const Mesh<T,NDim>& mesh ) : _mesh(mesh) {
      for( auto& comp : _comps ) {
        comp.resize( _mesh.size() );
        for ( auto& elm : comp ) elm = static_cast<T>(0);
        comp.shrink_to_fit();
      }
    }

    template < int Comp >
    Field& set_offset( const std::array< offset_t, DGrid >& offset ) {
      std::get<Comp>(_comps).offset = offset;
      return *this;
    }

    inline const auto& mesh() const { return _mesh;}

    template < int Comp, typename... Indices >
    inline auto& c( const Indices&... indices  ) {
      return std::get<Comp>(_comps)[_mesh.linearized_index_of_whole_mesh(indices)];
    }

    template < int Comp, typename... Indices >
    inline auto c( const Indices&... indices  ) {
      return std::get<Comp>(_comps)[_mesh.linearized_index_of_whole_mesh(indices)];
    }

  };
}

#endif
