#ifndef  INCLUDE_PREDEFS_H_
#define  INCLUDE_PREDEFS_H_

#include <math.h>
#include <limits>

#define VECTOR_DIM 3
// #define CONST_PI acos(-1.0)
#define CONST_PI 3.14159265358979323846
#define MAX_CELL std::numeric_limits<int>::max()
#define MAX_RC std::numeric_limits<float>::max()
#define MAX_UINT std::numeric_limits<unsigned int>::max()
#define EPS 1.0e-12
#define NUM_PTC_SPECIES 3
#define NUM_BOUNDARIES 6

#define CUDA_PTC_STREAMS 3

#define GUARD_TILES 27
#define REGIONS_PER_DIM 3
#define CENTRAL_REGION 14
#define MAX_PTC_PER_SPLIT 16384

typedef double Scalar;

#endif
