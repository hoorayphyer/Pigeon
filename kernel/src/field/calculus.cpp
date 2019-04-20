#include "field/calculus.hpp"

namespace field {
  // NOTE i is the linear index of the resultant field
  template < int Order, int IGrid, typename T, int DGrid >
  constexpr T D ( const Component<T,DGrid>& f, int i, int s ) noexcept {
    if constexpr ( IGrid >= DGrid ) return 0.0;
    const bool o = !std::get<IGrid>(f.offset); // offset of resultant field
    else if ( Order == 2 ) {
      return f[i + o * s ] - f[i - !o * s ];
        }
    static_assert( Order == 2 );
  }
}

namespace field {
  namespace impl {
    template < typename Grid, typename Extent >
    struct inbulk_t {
      const Grid& grid;
      const Extent& ext;

      template < int I >
      inline bool in_bulk( int i ) noexcept {
        return ( std::get<I>(grid).guard - 1 < i ) &&
          ( i < std::get<I>(ext) - std::get<I>(grid).guard );
      }
    };
  }


  // TODO this version has not treated boundary issues
  template < typename T, int DGrid, int Order >
  struct calc<knl::coordsys_t::Cartesian, T, int DGrid> {
    void curl( Field<T,3,DGrid>& output, const Field<T,3,DGrid>& input ) {
      const auto& stride = output.stride;
      auto f = inbulk_t{_grid, output.extent};
      int ijk[DGrid];
      for ( int I = 0; I < stride.back();  ) {
        ijk[0] = I % stride[1];
        if ( !f.in_bulk<0>(ijk[0]) ) { ++I; continue; }

        if constexpr ( DGrid >= 2 ) {
            if ( ijk[0] != grid[0].guard ) goto CURL;
            ijk[1] = ( I % stride[2] ) / stride[1];
            if ( !f.in_bulk<1>(ijk[1]) ) { I+=stride[1]; continue; }
          }

        if constexpr ( DGrid == 3 ) {
            if ( ijk[1] != grid[1].guard ) goto CURL;
            ijk[2] = I / stride[2];
            if ( !f.in_bulk<2>(ijk[2]) ) { I+=stride[2]; continue; }
          } static_assert( DGrid < 4 );

      CURL:
        output.c<0>[I] = D<DGrid>(1,field<2>) - drv<DGrid>(2,field<1>);


      }
    }

    void div( Field<T,1,DGrid>& output, const Field<T,3,DGrid>& input ) {
      
    }

    void lapl( Field<T,3,DGrid>& output, const Field<T,3,DGrid>& input );
  };








  void curl( const vector_field& input, vector_field& output,
             FieldType type, const bool isBoundary[],
             const Index& start, const Extent& ext ) {
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

  }
}
