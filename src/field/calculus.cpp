#include "field/calculus.hpp"

namespace field {
  template < typename Real, int Order >
  constexpr Real diff (Real f[], Real delta) noexcept {
    if constexpr ( Order == 2 ) return (f[1] - f[0]) / delta;
    static_assert( Order == 2 );
  }
}

namespace field {
  void curl( const vector_field& input, vector_field& output,
             FieldType type, const bool isBoundary[],
             const Index& start, const Extent& ext ) {
  auto time_start = high_resolution_clock::now();

  // Compute the 3 components of the result separately
  for (int i = 0; i < VECTOR_DIM; i++) {
    output.data(i).assign(0.0);
    // (Curl F)_i = D_j F_k - D_k F_j
    int j = (i + 1) % VECTOR_DIM;
    int k = (i + 2) % VECTOR_DIM;

    // First find the stagger of the input field type
    Index stagger_j = GetStagProperty(type, j);
    Index stagger_k = GetStagProperty(type, k);

    DiffParams mod;
    // Compute D_j F_k, dealing with singularity
    if (_dim > j) {
      mod.mod_bulk[0] = mod.mod_lower[0] = mod.mod_upper[0] = _h[k];
      mod.mod_bulk[1] = mod.mod_lower[1] = mod.mod_upper[1] = _h[j] | _h[k]; // We can also use +?
      // mod.assume_zero =
      if (j == 1 && k != 0 && (_coord_type == CoordType::SPHERICAL || _coord_type == CoordType::LOG_SPHERICAL)) {
        mod.mod_lower[0] = mod.mod_upper[0] = 0;
        mod.mod_lower[1] = mod.mod_upper[1] = _h[j];
        mod.factor[0] = mod.factor[2] = 2.0;
        mod.assume_zero[0] = mod.assume_zero[1] = true;
      }
      Derivative(input.data(k), output.data(i), j, stagger_k, start, ext,
                 isBoundary[2*j], isBoundary[2*j + 1], mod);
    }

    if (_dim > k) {
      // Compute D_k F_j, dealing with singularity
      mod.mod_bulk[0] = mod.mod_lower[0] = mod.mod_upper[0] = _h[j];
      mod.mod_bulk[1] = mod.mod_lower[1] = mod.mod_upper[1] = _h[j] | _h[k]; // We can also use +?
      mod.factor[0] = mod.factor[1] = mod.factor[2] = -1.0;
      if (k == 1 && j != 0 && (_coord_type == CoordType::SPHERICAL || _coord_type == CoordType::LOG_SPHERICAL)) {
        mod.mod_lower[0] = mod.mod_upper[0] = 0;
        mod.mod_lower[1] = mod.mod_upper[1] = _h[k];
        mod.factor[0] = mod.factor[2] = -2.0;
        mod.assume_zero[0] = mod.assume_zero[1] = true;
      }
      Derivative(input.data(j), output.data(i), k, stagger_j, start, ext,
                 isBoundary[2*k], isBoundary[2*k + 1], mod);
    }
  }

  auto time_end = high_resolution_clock::now();
}
}
