#ifndef INCLUDE_SCALES_H_
#define INCLUDE_SCALES_H_

#include "Grid.h"
#include "MultiArray.h"
#include "Predefs.h"

class Scales {
 private:
  // enum { stagger_num = (1 << DIM) } ;
  typedef MultiArray<Scalar> array_type;
  typedef std::vector<array_type> scale_array;

  std::array<scale_array, VECTOR_DIM> _scales;
  int _dim;
  int _stagger_num;

 public:
  Scales() : _scales() {}

  template <typename CoordSystem>
  Scales(const CoordSystem& coord, const Grid& grid) {
    _dim = grid.dimension;
    _stagger_num = (1 << _dim);
    for (int i = 0; i < VECTOR_DIM; i++) {
      _scales[i].resize(_stagger_num);
      for (int j = 0; j < _stagger_num; j++) {
        _scales[i][j].resize(grid.extent());
      }
    }

    for (int n = 0; n < _stagger_num; n++) {
      // Using bit operations might be more intuitive?
      int istag = check_bit(n, 0);  // Extract lowest bit
      int jstag = check_bit(n, 1);  // Extract second lowest bit
      int kstag = check_bit(n, 2);  // Extract highest bit
      // Alternatives
      // int istag = n % 2;
      // int jstag = (n / 2) % 2;
      // int kstag = (n / 4) % 2;

      // loop over cell indices(i,j,k)
      for (int k = 0; k < grid.dims[2]; k++) {
        for (int j = 0; j < grid.dims[1]; j++) {
          for (int i = 0; i < grid.dims[0]; i++) {
            // calculate the coordinate values
            double q1 = grid.pos(0, i, istag) + EPS;
            double q2 = grid.pos(1, j, jstag) + EPS;
            double q3 = grid.pos(2, k, kstag) + EPS;
            // calculate the scaling functions h1, h2, h3
            _scales[0][n](i, j, k) = coord.h1(q1, q2, q3);
            _scales[1][n](i, j, k) = coord.h2(q1, q2, q3);
            _scales[2][n](i, j, k) = coord.h3(q1, q2, q3);
          }
        }
      }
    }
  }

  ~Scales() {}

  array_type& operator()(int n, Vec3<int> stagger) {
    int stag_id = stagger.x;
    if (_dim >= 2) stag_id += stagger.y * 2;
    if (_dim >= 3) stag_id += stagger.z * 4;
    return _scales[n][stag_id];
  }

  const array_type& get_scales(int n, int stag_id) const {
    return _scales[n][stag_id];
  }

  array_type& get_scales(int n, int stag_id) { return _scales[n][stag_id]; }
};  // ----- end of class Scales -----

#endif  // INCLUDE_SCALES_H_
