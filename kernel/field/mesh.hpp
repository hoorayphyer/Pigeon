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

        constexpr bool operator!= ( const apt::BlockIteratorEnd& end ) const noexcept {
          return _bitr != end;
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
      // int I = i_bulk[NDim-1] + _margin[NDim-1][0] + _guard;
      int I = i_bulk[NDim-1] + _guard; // TODO check RVO conditions
      for ( int i = NDim - 2; i > -1; --i )
        // I = ( i_bulk[i] + _margin[i][0] + _guard ) + I * _extent[i];
        I = ( i_bulk[i] + _guard ) + I * _extent[i];
      return I;
    }

    // NOTE I_bulk_begin[ith_dim] is not significant
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

    constexpr const auto& extent() const noexcept { return _extent; }

    constexpr auto bulk_dim( int ith_dim ) const noexcept {
      // TODOL check bounds
      return _extent[ith_dim] - 2 * _guard;
    }

    constexpr apt::Index<D> bulk_dims() const noexcept {
      apt::Index<D> res;
      for ( int i = 0; i < D; ++i )
        res[i] = bulk_dim(i);
      return res;
    }

    // NOTE stride can be used to extract total linear size of the mesh by passing ith_dim = DGrid
    constexpr auto stride( int ith_dim ) const noexcept {
      // TODOL check bounds
      if ( 0 == ith_dim ) return 1;
      else return _extent[ith_dim - 1] * stride( ith_dim - 1 );
    }

    // // TODOL this potentially will introduce subtle bugs when using trI from the full mesh. Naturally the goal here is to have a mesh that is 1 dim short. But such a mesh doesn't work well with trI from the full mesh. One solution is to provide trI to int by subtracting g * stride_normal
    // auto squeeze( int ith_dim ) const noexcept {
    //   // TODOL check bounds
    //   apt::Index<NDim> bulk;
    //   for ( int i = 0; i < NDim; ++i ) bulk[i] = bulk_dim(i);
    //   bulk[ith_dim] = 1;
    //   // NOTE guard is still kept in the squeezed direction to conserve transindex compatibility
    //   return Mesh( bulk, guard() );
    // }

  };
}

#endif
