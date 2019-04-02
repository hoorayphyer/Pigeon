#ifndef  _FIELD_MESH_HPP_
#define  _FIELD_MESH_HPP_

#include "apt/block.hpp"
#include "apt/pair.hpp"

namespace field {
  template < int D >
  struct Mesh {
  private:
    apt::Index<D> _extent{}; // the full mesh, consisting of bulk extent and guard

    // TODOL define _margin for applying boundary conditions. May need to modify linearized_index and origin
    // apt::array< apt::pair<int>, D > _margin {}; // margin is where boundary conditions are applied. They don't include guard cells
    int _guard = 0; // uniform value in all directions

  public:
    static constexpr int NDim = D;
    constexpr Mesh() = default;

    constexpr Mesh( apt::Index<D> bulk_extent, int guard ) noexcept
      : _extent(bulk_extent), _guard(guard) {
      apt::foreach<0,D>
        ( [&guard]( auto& ext ) { ext += 2 * guard; }, _extent );
    }

    struct ProjBlock {
    private:
      const Mesh<D>& _mesh;
      int _stride_normal_in_mesh;
      apt::Index<D> _transIb_bulk;
      apt::Block<D> _block;

      constexpr ProjBlock( const Mesh<D>& mesh, int ith_dim, apt::Index<D> I_bulk_begin, apt::Index<D> extent )
        : _mesh(mesh),
          _stride_normal_in_mesh(mesh.stride(ith_dim)),
          _transIb_bulk( [&I_bulk_begin, ith_dim](){ I_bulk_begin[ith_dim] = 0; return std::move(I_bulk_begin);}() ),
          _block( [&extent, ith_dim](){ extent[ith_dim] = 1; return std::move(extent); }() ) {}
    public:
      friend class Mesh<D>;

      struct TransIndex {
        const int transindex;
        const int stride_normal;

        constexpr int operator|( int i_bulk_normal ) const noexcept {
          return i_bulk_normal * stride_normal  + transindex;
        }
      };

      struct ProjBlockIterator {
      private:
        apt::BlockIterator<D> _bitr;
        const ProjBlock& _pb;

      public:
        using difference_type = void;
        using value_type = void;
        using reference = TransIndex;
        using pointer = void;
        using iterator_category = std::forward_iterator_tag;

        constexpr ProjBlockIterator( apt::BlockIterator<D> bitr, const ProjBlock& projblock ) noexcept
          : _bitr(std::move(bitr)), _pb( projblock ) {}

        constexpr bool operator!= ( const apt::Index<D>& idx ) const noexcept {
          return _bitr != idx;
        }

        constexpr ProjBlockIterator& operator++() noexcept {
          ++_bitr;
          return *this;
        }

        constexpr ProjBlockIterator operator++(int) noexcept {
          auto res = *this;
          ++(*this);
          return res;
        }

        constexpr reference operator*() noexcept {
          // NOTE linearized_index automatically absorbs guard * stride_normal into transindex
          return { _pb._mesh.linearized_index_of_whole_mesh(*( _bitr ) + _pb._transIb_bulk), _pb._stride_normal_in_mesh };
        }
      };

      friend class ProjBlockIterator;

      constexpr auto begin() const noexcept {
        return ProjBlockIterator( _block.begin(), *this );
      };

      constexpr auto end() const noexcept {
        return _block.end();
      }
    };

    constexpr int linearized_index_of_whole_mesh( const apt::Index<NDim>& i_bulk ) const noexcept {
      // TODO check bounds on i_bulk???
      int i = NDim - 1;
      // int I = i_bulk[i] + _margin[i][0] + _guard;
      int I = i_bulk[i] + _guard;
      for ( --i; i > -1; --i )
        // I = ( i_bulk[i] + _margin[i][0] + _guard ) + I * _extent[i];
        I = ( i_bulk[i] + _guard ) + I * _extent[i];
      return I;
    }

    constexpr auto project( int ith_dim, apt::Index<D> I_bulk_begin, apt::Index<D> extent ) const noexcept {
      return ProjBlock( *this, std::move(ith_dim), std::move(I_bulk_begin), std::move(extent) );;
    }


    // constexpr const auto& margin() const noexcept { return _margin; }
    constexpr int guard() const noexcept { return _guard; }

    constexpr apt::Index<NDim> origin() const noexcept { // the first mesh cell expressed in bulk_indices
      apt::Index<NDim> origin;
      for ( int i = 0; i < NDim; i++ )
        // origin[i] = - _margin[i][0] - _guard;
        origin[i] = - _guard;
      return origin;
    }

    // constexpr int extent( int ith_dim ) const noexcept { // full extent of the mesh
    //   return _bulk[ith_dim].dim() + _margin[ith_dim][0] + _margin[ith_dim][1] + 2 * _guard;
    // }

    constexpr const auto& extent() const noexcept { return _extent; }

    constexpr auto bulk_dim( int ith_dim ) const noexcept {
      // TODOL check bounds
      return _extent[ith_dim] - 2 * _guard;
    }

    constexpr auto stride( int ith_dim ) const noexcept {
      // TODOL check bounds
      if ( 0 == ith_dim ) return 1;
      else return _extent[ith_dim - 1] * stride( ith_dim - 1 );
    }

  };
}

#endif
