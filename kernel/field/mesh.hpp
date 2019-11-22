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
    constexpr auto linear_size() const noexcept { return _stride.back(); }
    constexpr const auto& range() const noexcept { return _range; }
    constexpr const auto& range(int i) const noexcept { return _range[i]; }

    constexpr int linear_index( const apt::Index<D>& I ) const noexcept {
      // TODO check that _range.far_begin <= I < _range.far_end. NOTE that I failing to to so may still result in a valid linear index;
      // for ( int i = 0; i < D; ++i ) {
      //   if ( I[i] < _range[i].far_begin() or I[i] >= _range[i].far_end() ) {
      //     std::cout << "  I = " << I[0] << ", " << I[1] << ", range = " << _range[0].far_begin() << ", " << _range[0].far_end() << ", " << _range[1].far_begin() << ", " << _range[1].far_end() << std::endl;
      //     throw std::runtime_error("linear_index bad");
      //   }
      // }
      int res = _linear_offset;
      for ( int i = 0; i < D; ++i ) res += I[i] * _stride[i];
      return res;
    }

    constexpr int linear_index( const apt::Longidx& i ) const noexcept {
      // TODO check bounds on i.dir()
      return _linear_offset + i.val() * _stride[i.dir()];
    }

  };
}

#endif
