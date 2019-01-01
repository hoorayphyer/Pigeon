#ifndef _GRID_PARTICLE_INTERPLAY_HPP_
#define _GRID_PARTICLE_INTERPLAY_HPP_

namespace gpi {
  void interpolate_field ( Vec3 q, const Field& field, const Grid& grid, const ShapeFunction& sf ) {
    q -= field.shift;
    q /= grid.deltas;
    int ilb = grid.index(q - radius) + 1;
    int iub = grid.index(q + radius);
    Vec3 result;
    for ( ilb.k - ub.k )
      for ( ilb.j - ub.j )
        for ( ilb.i - ub.i )
          result += field(i,j,k) sf( di, dj, dk );
    return result
  }

}


#endif
