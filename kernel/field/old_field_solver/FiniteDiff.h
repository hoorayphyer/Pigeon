#ifndef INCLUDE_FINITEDIFF_H_
#define INCLUDE_FINITEDIFF_H_

#include <array>
#include <cstdint>
#include <vector>

#include "Fields.h"

//#define ORDER 4

enum class CoordType;
class FieldCommunicator;
class Scales;

struct DiffParams {
  char mod_lower[2] = {0, 0};  // { mult_mod, div_mod }
  char mod_bulk[2] = {0, 0};
  char mod_upper[2] = {0, 0};
  Scalar factor[3] = {1.0, 1.0, 1.0};
  bool assume_zero[2] = {false, false};
};

class FiniteDiff {
 private:
  typedef MultiArray<Scalar> array_type;
  //  typedef FiniteDiff<Coord> self_type;

  const CoordType _coord_type;

  Grid _grid;
  // intermediate fields used in ComputeLaplacian
  ScalarField<Scalar> _scalar_tmp;
  VectorField<Scalar> _vector_tmp;
  Scales* _pScales;
  int _dim;
  int _order;
  int _h[3] = {1, 2, 4};
  // Timer _timer;

  // default computation domain
  Index _d_start0;
  Extent _d_ext0;

  void Derivative(const MultiArray<Scalar>& input, MultiArray<Scalar>& output,
                  int dir, const Index& stagger, bool lower, bool upper,
                  const DiffParams& mod);

  void Derivative(const MultiArray<Scalar>& input, MultiArray<Scalar>& output,
                  int dir, const Index& stagger, const Index& start,
                  const Extent& ext, bool lower, bool upper,
                  const DiffParams& mod);

  // void GetDefExt(Index& start, Extent& ext, const bool isBoundary[]);

 public:
  typedef VectorField<Scalar> vector_field;
  typedef ScalarField<Scalar> scalar_field;

  FiniteDiff(CoordType coordtype, const Grid& grid, int order = 2);
  virtual ~FiniteDiff();

  void SetDefaultComputeDomain(Index d_start, Extent d_ext);
  // self_type& operator=(const self_type&);
  // self_type& operator=(self_type&&);
  void ComputeDivergence(const vector_field& input, scalar_field& output,
                         FieldType input_type, const bool isBoundary[]);
  void ComputeDivergence(const vector_field& input, scalar_field& output,
                         FieldType input_type, const bool isBoundary[],
                         const Index& start, const Extent& ext);

  void ComputeGradient(const scalar_field& input, vector_field& output,
                       StaggerType input_type, const bool isBoundary[]);
  void ComputeGradient(const scalar_field& input, vector_field& output,
                       StaggerType input_type, const bool isBoundary[],
                       const Index& start, const Extent& ext);

  void ComputeCurl(const vector_field& input, vector_field& output,
                   FieldType input_type, const bool isBoundary[]);
  void ComputeCurl(const vector_field& input, vector_field& output,
                   FieldType input_type, const bool isBoundary[],
                   const Index& start, const Extent& ext);

  void ComputeLaplacian(const vector_field& input, vector_field& output,
                        FieldType input_type, const bool isBoundary[],
                        FieldCommunicator& fc, bool skipDiv = false);

  void ComputeLaplacian(const vector_field& input, vector_field& output,
                        FieldType input_type, const bool isBoundary[],
                        FieldCommunicator& fc, const Index& start,
                        const Extent& ext, bool skipDiv = false);

  array_type& get_scales(int n, int istag, int jstag, int kstag);

  const Scales* scales() const { return _pScales; }
};  // ----- end of class FiniteDiff -----

#endif  // INCLUDE_FINITEDIFF_H_
