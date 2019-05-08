#ifndef _APERTURE_PARTICLE_BC_HPP_
#define _APERTURE_PARTICLE_BC_HPP_

#include "field/field.hpp"
#include <cstdlib> // std::labs

namespace aperture {
  // NOTE by design WJ has offset = MIDWAY in all directions
  template < typename T, int DGrid >
  struct PBC_WJ_FoldBackAtAxis {
    // TODO need calibrator between gI_start and field. This can get rid of the margin in mesh
    // both gI_start_incl and gI_end_incl are inclusive. This is to maintain the symmetry between upper and lower bounds as much as possible
    void apply( field::Field<T,3,DGrid>& WJ, int normal, const apt::Index<DGrid> gI_start_incl, apt::Index<DGrid> gI_end_incl ) {
      bool is_lb = gI_start[normal] <= gI_end[normal]; // check if at lower bound. Only the order in the normal direction matters
      auto& extent = gI_end_incl; // recycle
      for ( int i = 0; i < DGrid; ++i )
        extent[i] = std::labs( gI_end_incl[i] - gI_start_incl[i] );

      int ext_normal = extent[normal];
      extent[normal] = 1; // make extent cover only the transverse directions
      for ( int i = 0; i < ext_normal; ++i ) {
        for ( auto I : apt::block(extent) ) {
          I[normal] = i;
          I += gI_start;
          auto I_mirror = I;
          I_mirror[normal] = 2 * gI_start[normal] + ( is_lb ? -1 : 1) - I[normal];
          WJ[0](I) += WJ[0](I_mirror);
          WJ[1](I) += WJ[1](I_mirror); // TODO check this!! old code has +=, this conserves current when doing the scan. NOTE WJ[i where i < DGrid] is not current yet
          WJ[2](I) -= WJ[2](I_mirror);

          WJ[0](I_mirror) = 0.0;
          WJ[1](I_mirror) = 0.0;
          WJ[2](I_mirror) = 0.0;
        }
      }
    }

  };

}

#endif
