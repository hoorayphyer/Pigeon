#ifndef  _FIELD_MESH_HPP_
#define  _FIELD_MESH_HPP_

#include "apt/block.hpp"
#include "apt/range.hpp"

// NOTE all indices in the function call interfaces are domain indices
namespace field {
  template < int D >
  struct Mesh {
  private:
    apt::array<apt::Range,D> _range;
    apt::array<int,D+1> _stride;
    int _linear_offset = 0;

  public:
    static constexpr int NDim = D;
    constexpr Mesh() = default;

    constexpr Mesh( const apt::array<apt::Range,D>& range ) noexcept
      : _range(range) {
      _stride[0] = 1;
      for ( int i = 0; i < D; ++i ) _stride[i+1] = _stride[i] * (_range[i].full_size());
      for ( int i = 0; i < D; ++i ) _linear_offset -= _stride[i] * ( _range[i].far_begin());
    }

    constexpr const auto& stride() const noexcept { return _stride; }
    constexpr const auto& range() const noexcept { return _range; }
    constexpr const auto& range(int i) const noexcept { return _range[i]; }

    constexpr int linear_index( const apt::Index<D>& I ) const noexcept {
      // TODO check bounds on I
      int res = _linear_offset;
      for ( int i = 0; i < D; ++i ) res += I[i] * _stride[i];
      return res;
    }

    constexpr int linear_index( const apt::Longidx& i ) const noexcept {
      // TODO check bounds on I
      return _linear_offset + i.val() * _stride[i.dir()];
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
