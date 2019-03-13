#include "FiniteDiff.h"
#include "CoordSystem.h"
#include "ArrayOperations.h"
#include "Scales.h"
#include "Domain.h"

#define H1 1
#define H2 2
#define H3 4

inline int flip (int n) {
  return (n == 1 ? 0 : 1);
}

inline void multiplyMod (char mod_mult, Scalar& value, Scales* scales, int idx, const Index& stagger) {
  if ((mod_mult & H1) == H1) value *= (*scales)(0, stagger)[idx];
  if ((mod_mult & H2) == H2) value *= (*scales)(1, stagger)[idx];
  if ((mod_mult & H3) == H3) value *= (*scales)(2, stagger)[idx];
}

inline void divideMod (char mod_div, Scalar& value, Scales* scales, int idx, const Index& stagger) {
  if ((mod_div & H1) == H1) value /= (*scales)(0, stagger)[idx];
  if ((mod_div & H2) == H2) value /= (*scales)(1, stagger)[idx];
  if ((mod_div & H3) == H3) value /= (*scales)(2, stagger)[idx];
}

inline Scalar diff_2nd (Scalar f[], Scalar delta) {
  return (f[1] - f[0]) / delta;
}

inline Scalar diff_2nd_side (Scalar f[], Scalar delta, bool zero) {
  return (zero ? (3.0 * f[0] - f[1] / 3.0) / delta :
          (-2.0 * f[0] + 3.0 * f[1] - f[2]) / delta);
}

// inline Scalar diff_2nd_side (Scalar f[], Scalar delta) {
//   return (3.0 * f[0] - f[1] / 3.0) / delta;
// }

inline Scalar diff_4th (Scalar f[], Scalar delta) {
  return ((f[2] - f[1]) * 9.0 / 8.0 - (f[3] - f[0]) / 24.0) / delta;
}

// this is the (1, rest)-type derivative
inline Scalar diff_4th_side (Scalar f[], Scalar delta) {
  return (-f[0] * 11.0 / 12.0 + f[1] * 17.0 / 24.0 + f[2] * 3.0 / 8.0 - f[3] * 5.0 / 24.0 + f[4] / 24.0) / delta;
}

// this is the (0, all)-type derivative.
inline Scalar diff_4th_side_mod (Scalar f[], Scalar delta, bool zero) {
  return (zero ? (f[0] * 35.0 / 8.0 - f[1] * 35.0 / 24.0 + f[2] * 21.0 / 40.0 - f[3] * 5.0 / 56.0) / delta :
          (-f[0] * 31.0 / 8.0 + f[1] * 229.0 / 24.0 - f[2] * 75.0 / 8.0 + f[3] * 37.0 / 8.0 - f[4] * 11.0 / 12.0) / delta);
}

// inline Scalar diff_4th_side_mod_zero (Scalar f[], Scalar delta) {
//   return (f[0] * 35.0 / 8.0 - f[1] * 35.0 / 24.0 + f[2] * 21.0 / 40.0 - f[3] * 5.0 / 56.0) / delta;
// }

/// Apply derivative to each point of the grid using a loop
void
FiniteDiff::Derivative(const MultiArray<Scalar>& input, MultiArray<Scalar>& output,
                             int dir, const Index& stagger, const Index& start, const Extent& ext,
                             bool lower, bool upper, const DiffParams& mod) {
  // Explanation of Derivative
  // 0. Derivative calculates the derivative of input and ADD it to output.
  // 1. the passed in start and ext are the computational domain for an input field unstaggered in all directions. The staggeredness is taken care of in the Derivative automatically.
  // 2. the computational domain being defined as above, we always requires that all fields that lie on the edges of this domain also belong to the domain.
  // 3. range_start and range_end are adjusted for the output field. Hence so is the corresponding loop index i.
  // 4. bdry_layer is used to specify the number of cells of output field near the boundary for which some kind of one_sided derivative should be used. Note it bdry_layer depends on input stagger. bdry_layer is anyway calculated, but is only effective at true boundaries.
  // 5. mod_lower[2] = { mod_mult, mod_div }. These are relevant only in the case of evaluating an output field exactly on the boundary. So if the input field is staggered, or if an output field is evaluated at a cell in bdry_layer other than the immediate cell on the boundary, mod_bulk should be used.
  // 6. diff_2nd and diff_2nd_side all take in the address of the field with smallest index. At upper boundary though, for one_sided diff, it is the field with highest index that's passed in, so that the same diff_side function can be reused with delta being negative. Same applies to 4th order.
  // 7. cached values are of length blocksize + 2 * guard. It is important to note which is the first value cached given any result range ( i.e. the range in which results are evaluated ). For bulk, the first value cached is the leftmost value in the "guard cells" of the result range. In the code, i - _order/2 is used to refer to this cell. For lower and upper boundaries, it is always the first cell beyond the boundary that is cached, regardless of guard. The need to include one extra cell is because in higher than 2nd order diff schemes at the lower boundary, one-sided diff for staggered input fields will need the value right on the boundary. For upper boundary, the first cached value is actually unnecessary; the stagger we work with guarantees that at upper boundary, no values beyond the boundary is needed. But this is the convention.
  // 8. there is a variable called range when caching values on boundaries. It is equal to the largest possible number of input fields that is needed to evaluate one-sided derivative. Note this number doesn't any input field beyond the boundary. In other words, this range is equal to the number of input fields when result is evaluated right on the boundary. I don't know why it is range+2 number of fields points cached; seems that range+1 or even just range is sufficient.
  // 9.[Problem] cached values at boundary seems a bit strange. The implementation actually keeps caching for each one of the boundary cell, which seems unnecessary. This is not a problem for 2nd order.
  // 10. range_start in transverse directions is modified based on the staggeredness of field in that direction, no matter whether it's at true lower boundary or not. This is the right behavior if at the true boundary in that transverse direction. In bulks, this value will be equal to and eventually be overwritten by the sendGuardCells values.

  if (dir < 0 || dir > 2)
    throw std::invalid_argument("Direction has to be 0 to 2!");

  const int block_size = 8;

  // First try to obtain the range of the loop from given start and
  // extent parameters
  Vec3<int> range_start { std::max(_grid.guard[0], start[0]),
        std::max(_grid.guard[1], start[1]),
        std::max(_grid.guard[2], start[2]) };
  Vec3<int> range_end { std::min(_grid.guard[0] + _grid.reducedDim(0), start[0] + ext[0]),
        std::min(_grid.guard[1] + _grid.reducedDim(1), start[1] + ext[1]),
        std::min(_grid.guard[2] + _grid.reducedDim(2), start[2] + ext[2]) };
  // Modify the range depending on stagger. On directions
  // perpendicular to derivative direction, we include an extra cell
  // if the stagger in that direction is 1. On direction of the
  // derivative, we do the opposite because the result has different
  // stagger from the input stagger. We also obtain the transverse
  // directions on the side.
  int trans[2];
  int c = 0;
  for (int i = 2; i >= 0; i--) {
    if (dir != i) {
      trans[c] = i;
      c += 1;
      if (i < _dim) // Don't want to mess with untracked dimensions
        range_start[i] -= stagger[i];
    } else {
      if (lower)
        range_start[i] -= flip(stagger[i]);
    }
  }

  // Determine the thickness of boundary layer
  // When stagger[dir] = 1 and at least one end assumes zero,
  // it still needs special treatment.
  // FIXME: maintain bndry_layers separately for lower and upper?
  int bndy_layer = 0;
  if (_order == 2) bndy_layer = 1 - stagger[dir];
//          * static_cast<int>( !( mod.assume_zero[0] ) && !( mod.assume_zero[1] ) );
  else if (_order == 4) bndy_layer = 2 - stagger[dir];
//          * static_cast<int>( !( mod.assume_zero[0] ) && !( mod.assume_zero[1] ) );

  // We need the stagger configuration of the result to compute
  // mod_div
  Index stagger_result = stagger;
  stagger_result[dir] = flip(stagger[dir]);

  // This is the buffer array for one derivative in a block

  // Outer loop goes over transverse directions. Note: The indices i,
  // j, k correspond to the indices of the result field point.
// #pragma omp parallel
  // {
    // for collapse(2)
    Scalar* values = new Scalar[block_size + 2 * _grid.guard[dir]];
    // if (omp_get_thread_num() == 0)
    //   std::cout << "Openmp thread number is " << omp_get_num_threads() << std::endl;

// #pragma omp for collapse (2)
    for (int k = range_start[trans[0]]; k < range_end[trans[0]]; k++) {
      for (int j = range_start[trans[1]]; j < range_end[trans[1]]; j++) {
        // Can't use loop index to deal with boundary in transverse
        // direction, need to check for singularities

        // This is the transverse index
        int transIdx = j * DirIncrement(trans[1], _grid.extent()) + k * DirIncrement(trans[0], _grid.extent());

        // Inner loop goes over the derivative direction. Take
        // block_size strides to minimize memory operations.

        int i = range_start[dir];
        char mod_mult = mod.mod_lower[0], mod_div = mod.mod_lower[1], factor = mod.factor[0];
        ///////////////////////////////////////////////////////////////////////
        // Lower boundary. Since i was defined and initialized above,
        // skip initialization
        ///////////////////////////////////////////////////////////////////////
        for ( ; i < range_start[dir] + (bndy_layer * (int)lower); i++) {
          if (i > range_start[dir] || stagger[dir] != 0) {
            mod_mult = mod.mod_bulk[0];
            mod_div = mod.mod_bulk[1];
            factor = mod.factor[1];
          }

          int range = 0;
          if (_order == 2) range = 3;
          else if (_order == 4) range = 5;

          for (int n = 0; n < range + 2; n++) {
            int idx = (i + n - stagger[dir]) * DirIncrement(dir, _grid.extent()) + transIdx;
            values[n] = input[idx];
            multiplyMod(mod_mult, values[n], _pScales, idx, stagger);
          }

          // adjust values[n] for staggered input field when zero is assumed.
//          if( i == range_start[dir] && stagger[dir] == 1 && mod.assume_zero[0] )
//              values[0] = 0.0;

          // if (i == 2)

          Scalar result = 0.0;
          // This is the index for the output field
          int outIdx = i * DirIncrement(dir, _grid.extent()) + transIdx;
          if (_order == 2) {
//            if ( stagger[dir] == 0) {
              result = factor * diff_2nd_side (values + 1 - stagger[dir], _grid.delta[dir], mod.assume_zero[0]);
//            } else {
//              result = factor * diff_2nd( values + 1 -stagger[dir], _grid.delta[dir]);
//            }
          } else if (_order == 4) {
            if (i == range_start[dir] && stagger[dir] == 0) {
              result = factor * diff_4th_side_mod (values + 1 - stagger[dir], _grid.delta[dir], mod.assume_zero[0]);
            } else {
              result = factor * diff_4th_side (values, _grid.delta[dir]);
            }
          }
          divideMod(mod_div, result, _pScales, outIdx, stagger_result);
          output[outIdx] += result;
        }

        int bulk_end = range_end[dir] - (bndy_layer * (int)upper);
        // Set modifier to be bulk
        mod_mult = mod.mod_bulk[0];
        mod_div = mod.mod_bulk[1];
        factor = mod.factor[1];
        ///////////////////////////////////////////////////////////////////////
        // Bulk
        // Again skip initialization because we reuse the same index
        ///////////////////////////////////////////////////////////////////////
        for ( ; i < bulk_end ; i += block_size) {
          // Fill the buffer array with values first
          for (int n = 0; n < block_size + 2 * _grid.guard[dir]; n++ ) {
            // if (i - _grid.g)
            // FIXME: _order / 2 is a good indicator?
            if (i - _order/2 + n > _grid.dims[dir] - 1) continue;
            // printf ("Exceeding limit: i is %d, n is %d!\n", i, n);
            int idx = std::min(i - _order/2 + n, _grid.dims[dir] - 1) * DirIncrement(dir, _grid.extent()) + transIdx;
            // if (idx >= input.size())
            //   Logger::err("Exceeding boundary at i =", i, "n = ", n);
            // else if (idx < 0)
            //   Logger::err("index is smaller than zero at i =", i, "n = ", n);
            values[n] = input[idx];
            multiplyMod(mod_mult, values[n], _pScales, idx, stagger);
          }
          for (int n = 0; n < block_size; n++) {
            if (i + n >= bulk_end) break;
            // This is the index for the output field
            int outIdx = (i + n) * DirIncrement(dir, _grid.extent()) + transIdx;
            // if (outIdx >= _grid.size())
            //   printf("!!!Grid size exceeded at i = %d, n = %d\n", i, n);
            // output[outIdx] = 0.0;
            Scalar result = 0.0;
            if (_order == 2)
              result = factor * diff_2nd (values + n + 1 - stagger[dir], _grid.delta[dir]);
            else if (_order == 4)
              result = factor * diff_4th (values + n + 1 - stagger[dir], _grid.delta[dir]);
            divideMod(mod_div, result, _pScales, outIdx, stagger_result);
            output[outIdx] += result;
          }
        }

        i = bulk_end;

        ///////////////////////////////////////////////////////////////////////
        // Upper boundary
        ///////////////////////////////////////////////////////////////////////
        for (; i < range_end[dir]; i++) {
          if (i == range_end[dir] - 1 && stagger[dir] == 0) {
            // Set modifiers to the upper boundary one
            mod_mult = mod.mod_upper[0];
            mod_div = mod.mod_upper[1];
            factor = mod.factor[2];
          }
          int range = 0;
          if (_order == 2) range = 3;
          else if (_order == 4) range = 5;

          for (int n = 0; n < range + 2; n++) {
            int idx = (i + 1 - stagger[dir] - n) * DirIncrement(dir, _grid.extent()) + transIdx;
            values[n] = input[idx];
            multiplyMod(mod_mult, values[n], _pScales, idx, stagger);
          }

          // adjust values[n] for staggered input field when zero is assumed.
//          if( i == range_end[dir] - 1 && stagger[dir] == 1 && mod.assume_zero[1] )
//              values[0] = 0.0;

          Scalar result = 0.0;
          // This is the index for the output field
          int outIdx = i * DirIncrement(dir, _grid.extent()) + transIdx;
          if (_order == 2) {
//            if ( stagger[dir] == 0 ) {
              result = factor * diff_2nd_side (values + 1 - stagger[dir], -_grid.delta[dir], mod.assume_zero[1]);
//            } else {
//              result = factor * diff_2nd( values + 1 - stagger[dir], -_grid.delta[dir]);
//            }
          } else if (_order == 4) {
            if (i == range_end[dir] - 1 && stagger[dir] == 0)
              result = factor * diff_4th_side_mod (values + 1 - stagger[dir], -_grid.delta[dir], mod.assume_zero[1]);
            else
              result = factor * diff_4th_side (values, -_grid.delta[dir]);
          }
          divideMod(mod_div, result, _pScales, outIdx, stagger_result);
          output[outIdx] += result;
        }
      }
    }

    delete[] values;
  // }
}

void
FiniteDiff::Derivative(const MultiArray<Scalar>& input, MultiArray<Scalar>& output, int dir, const Index& stagger,
                             bool lower, bool upper, const DiffParams& mod) {
  Derivative(input, output, dir, stagger, _d_start0, _d_ext0, lower, upper, mod);
}

FiniteDiff::FiniteDiff( CoordType coordtype, const Grid& grid, int order)
  : _coord_type(coordtype), _grid(grid),
    _scalar_tmp(grid), _vector_tmp(grid), _order(order) {

  _pScales = new Scales( ScalesLogSpherical(), grid);
  _dim = _grid.dimension;
  // initialize default computation domain
  _d_start0 = Index(_grid.guard[0], _grid.guard[1], _grid.guard[2]);
  _d_ext0 = Extent(_grid.reducedDim(0), _grid.reducedDim(1), _grid.reducedDim(2));

}

FiniteDiff::~FiniteDiff() {
  if (_pScales != nullptr) delete _pScales;
}

MultiArray<Scalar>&
FiniteDiff::get_scales(int n, int istag, int jstag, int kstag) {
  return (*_pScales)(n, Vec3<int>(istag, jstag, kstag));
}

void FiniteDiff::SetDefaultComputeDomain(Index d_start, Extent d_ext){
  _d_start0 = d_start;
  _d_ext0 = d_ext;
}

void
FiniteDiff::ComputeCurl(const vector_field& input, vector_field& output,
                        FieldType type, const bool isBoundary[],
                        const Index& start, const Extent& ext) {

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
      if (j == 1 && k != 0 && ( _coord_type == CoordType::LOG_SPHERICAL)) {
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
      if (k == 1 && j != 0 && ( _coord_type == CoordType::LOG_SPHERICAL)) {
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

void
FiniteDiff::ComputeCurl(const vector_field& input, vector_field& output,
                               FieldType type, const bool isBoundary[]) {
  ComputeCurl(input, output, type, isBoundary, _d_start0, _d_ext0);
}

void
FiniteDiff::ComputeDivergence(const vector_field &input, scalar_field &output,
                                     FieldType input_type, const bool *isBoundary,
                                     const Index& start, const Extent& ext) {
  output.assign(0.0);
  // Only loop over the first two dimensions if _dim is 2
  for (int i = 0; i < _dim; i++) {
    Index stagger = GetStagProperty(input_type, i);

    int j = (i + 1) % VECTOR_DIM;
    int k = (i + 2) % VECTOR_DIM;

    DiffParams mod;

    // Compute D_i F_i, dealing with singularity
    mod.mod_bulk[0] = mod.mod_lower[0] = mod.mod_upper[0] = _h[j] | _h[k];
    mod.mod_bulk[1] = mod.mod_lower[1] = mod.mod_upper[1] = _h[i] | _h[j] | _h[k];

    if (i == 1 && ( _coord_type == CoordType::LOG_SPHERICAL)) {
      mod.mod_lower[0] = mod.mod_upper[0] = 0;
      mod.mod_lower[1] = mod.mod_upper[1] = _h[i];
      mod.factor[0] = mod.factor[2] = 2.0;
      mod.assume_zero[0] = mod.assume_zero[1] = true;
    }
    Derivative(input.data(i), output.data(), i, stagger, start, ext,
               isBoundary[2*i], isBoundary[2*i + 1], mod);
  }
}

void
FiniteDiff::ComputeDivergence(const vector_field &input, scalar_field &output,
                                     FieldType input_type, const bool *isBoundary) {
  ComputeDivergence(input, output, input_type, isBoundary, _d_start0, _d_ext0);
}

void
FiniteDiff::ComputeGradient(const scalar_field &input, vector_field &output,
                                   StaggerType input_type, const bool *isBoundary,
                                   const Index& start, const Extent& ext) {
  // The gradient function always operates on a scalar field,
  // therefore it is a misnomer to specify input staggering as either
  // electric or magnetic field. It is better to specify it is
  // staggered or not, in all directions.
  int type = static_cast<int>(input_type);
  Index stagger = Index(type, type, type);
  output.assign(0.0);

  // Only loop over the first two dimensions if _dim is 2
  for (int i = 0; i < _dim; i++) {
    DiffParams mod;

    mod.mod_bulk[0] = mod.mod_lower[0] = mod.mod_upper[0] = 0;
    mod.mod_bulk[1] = mod.mod_lower[1] = mod.mod_upper[1] = _h[i];

//    if (i == 1 && (Coord == CoordType::SPHERICAL || Coord == CoordType::LOG_SPHERICAL)) {
//      mod.assume_zero[0] = mod.assume_zero[1] = true;
//    }
    // TODO: in 3D, deal with 1/h_phi at theta boundaries??

    Derivative(input.data(), output.data(i), i, stagger, start, ext,
               isBoundary[2*i], isBoundary[2*i + 1], mod);
  }
}

void
FiniteDiff::ComputeGradient(const scalar_field &input, vector_field &output,
                                   StaggerType input_type, const bool *isBoundary) {
  ComputeGradient(input, output, input_type, isBoundary, _d_start0, _d_ext0);
}

void
FiniteDiff::ComputeLaplacian(const vector_field &input, vector_field &output, FieldType input_type, const bool isBoundary[], Domain &domain, const Index& start, const Extent& ext, bool skipDiv) {
  // Obtain the correct stagger
  FieldType antitype = (FieldType::ETYPE == input_type)?
      FieldType::BTYPE : FieldType::ETYPE;
  StaggerType gradtype = (FieldType::ETYPE == input_type)?
      StaggerType::STAGGER_MINUS : StaggerType::STAGGER_PLUS;

  // Invoke the series of differential operators
  ComputeCurl(input,_vector_tmp,input_type,isBoundary, start, ext);
  domain.SendGuardCells(_vector_tmp);

  // Skip the divergence part if the flag is false
  if (!skipDiv) {
    ComputeDivergence(input,_scalar_tmp,input_type, isBoundary, start, ext);
    domain.SendGuardCells(_scalar_tmp);
  }

  // Adjust for conductor boundary where a shift of the evaluation
  // boundary of the second curl is required
  // if (pBC != nullptr) {
  //   for (int i = 0; i < NUM_BOUNDARIES; i++) {
  //     if (domain.isBoundary(i) && isConductor(pBC -> fieldBC((BoundaryPosition)i))) {
  //       if (i % 2 == 0) {
  //         start[i / 2] += 1;
  //         ext[i / 2] -= 1;
  //       } else {
  //         ext[i / 2] -= 1;
  //       }
  //     }
  //   }
  // }

  ComputeCurl(_vector_tmp,output,antitype,isBoundary, start, ext);

  if (!skipDiv) {
    ComputeGradient(_scalar_tmp,_vector_tmp,gradtype,isBoundary, start, ext);
  } else {
    _vector_tmp.assign(0.0);
  }

  // Put both components together
  for(int comp=0; comp < VECTOR_DIM; comp++){
    subtract(output.data(comp).begin(), _vector_tmp.data(comp).begin(), output.grid().extent());
    multiply(output.data(comp).begin(), output.grid().extent(), -1.0);
  }
}

void
FiniteDiff::ComputeLaplacian(const vector_field &input, vector_field &output, FieldType input_type, const bool isBoundary[], Domain &domain, bool skipDiv) {
  ComputeLaplacian(input, output, input_type, isBoundary, domain, _d_start0, _d_ext0, skipDiv);
}

// void
// FiniteDiff::ComputeLaplacian(const scalar_field &input, scalar_field &output, StaggerType input_type, Domain &domain) {}

// Instantiate the class with all the desired coordinate systems
//INSTANTIATE_CLASS_WITH_COORDTYPES(FiniteDiff);

