#ifndef FIELDBC_ROTATING_CONDUCTOR_H
#define FIELDBC_ROTATING_CONDUCTOR_H

#include "FieldBC.h"
#include "FUParams.h"

class FieldBC_rotating_conductor : public FieldBC {
private:
  BoundaryPosition _bpos;
  Grid _grid_cond;
  Index _shift; // shift stores the location of (0,0,0) of _grid_cond in the parent grid on the domain
  FBC::BCFunc_split_t _f_t; // time dependence of omega
  VectorField<Scalar> _E_cond;
  VectorField<Scalar> _B_cond;

  // for each field component, enforce the boundary condition in the boundary direction
  // from dir_start to dir_end
  void EnforceRotCond( Scalar ft, MultiArray<Scalar>& field_comp, const MultiArray<Scalar>& field_cond,
                       int dir_start, int dir_end ) const {
    int dir = static_cast<int>(_bpos) / 2;
    int trans[2] = { (dir + 1) % VECTOR_DIM , (dir + 2) % VECTOR_DIM };

    // since i, j, k are not necessarily 0, 1, 2, one needs to use an Index object for the actual indices of a field
    Index idx;

    for ( int k = 0; k < _grid_cond.dims[trans[1]]; ++k ) {
      idx[trans[1]] = k;
      for ( int j = 0; j < _grid_cond.dims[trans[0]]; ++j ) {
        idx[trans[0]] = j;
        for ( int i = dir_start; i < dir_end; ++i ) {
          idx[dir] = i;
          field_comp( idx[0] + _shift[0], idx[1] + _shift[1], idx[2] + _shift[2] )
              = ft * field_cond( idx[0], idx[1], idx[2] );
        }
      }
    }

  }

  void SetFieldCache( FieldType field_type, int comp, const FBC::BCFunc_split_x& fx ) {
    // although some values might not be used due to stagger
    // still calculate all values on _E_cond and _B_cond.
    Index stagger = GetStagProperty( field_type, comp );
    auto& field = FieldType::ETYPE == field_type ? _E_cond.data(comp) : _B_cond.data(comp);

    // note the case of being upper has been accommondated by _grid_cond.lower
    for ( int k = 0; k < _grid_cond.dims[2]; ++k ) {
      Scalar q3 = _grid_cond.pos( 2, k, stagger[2] );
      for ( int j = 0; j < _grid_cond.dims[1]; ++j ) {
        Scalar q2 = _grid_cond.pos( 1, j, stagger[1] );
        for ( int i = 0; i < _grid_cond.dims[0]; ++i ) {
          Scalar q1 = _grid_cond.pos( 0, i, stagger[0] );

          field(i,j,k) = fx(q1, q2, q3);
        }
      }
    }

  }

public:
  FieldBC_rotating_conductor( BoundaryPosition bpos, const Grid& grid, const FBC& fieldBC )
    : _bpos(bpos), _grid_cond(grid), _f_t(fieldBC.ft)  {
    int dir = bpos / 2;
    bool islower = bpos % 2;

    //1. taylor the full grid to fit the conductor region.
    //2. guard cells are irrelevant for _grid_cond. However, we still make
    // the grid_cond contain guard cells in the interested direction
    //3. the dimension of the grid will be indent. NOTE that the following applies to only boundaries in the q1 direction, or r in the case of log spherical; this depends on stagger. The principle is not to impose anything for components inside computational domain ( excluding edges ); on the edges, impose only if the components are continuous.
    _grid_cond.dims[dir] = grid.indent[bpos];
    _grid_cond.indent[dir * 2] = _grid_cond.guard[dir];
    _grid_cond.indent[dir * 2 + 1] = _grid_cond.guard[dir];

    // adjust grid start info based on lower or upper
    if( !islower ) {
      // note lower is defined to be at the right edge of the last guard cell.
      _grid_cond.lower[dir] = grid.pos( dir, grid.dims[dir] - _grid_cond.dims[dir], 1 );
      _shift[dir] = grid.dims[dir] - _grid_cond.dims[dir];
    }

    // resize reference fields
    _E_cond.resize( _grid_cond );
    _B_cond.resize( _grid_cond );

    // set up field caches
    SetFieldCache( FieldType::ETYPE, 0, fieldBC.E1 );
    SetFieldCache( FieldType::ETYPE, 1, fieldBC.E2 );
    SetFieldCache( FieldType::ETYPE, 2, fieldBC.E3 );

    SetFieldCache( FieldType::BTYPE, 0, fieldBC.B1 );
    SetFieldCache( FieldType::BTYPE, 1, fieldBC.B2 );
    SetFieldCache( FieldType::BTYPE, 2, fieldBC.B3 );
  }

  virtual void Apply(VectorField<Scalar>& Efield, VectorField<Scalar>& Bfield, Scalar time) override {
    // At lower boundary, discontious ( continuous ) components will be applied boundary conditions up to the indent-2-th ( indent-1-th ) cell.
    // At upper boundary, discontious and continuous components will both be applied boundary conditions starting from 0-th cell; stagger doesn't affect in this case.
    int dir = _bpos / 2;
    int trans[2] = { ( dir + 1 ) % VECTOR_DIM, ( dir + 2 ) % VECTOR_DIM };
    bool islower = _bpos % 2;

    Scalar ft = _f_t( time );

    // deal with continuous components. In this case, the full range of _grid_cond is used. Note the last parameter of EnforceRotCond is the index past actual end.
    EnforceRotCond( ft, Efield.data(trans[0]), _E_cond.data(trans[0]), 0, _grid_cond.dims[dir] );
    EnforceRotCond( ft, Efield.data(trans[1]), _E_cond.data(trans[1]), 0, _grid_cond.dims[dir] );
    EnforceRotCond( 1.0, Bfield.data(dir), _B_cond.data(dir), 0, _grid_cond.dims[dir] );

    // deal with discontinuous components. In this case, see comments in the beginning for adjustment details.Again, note the last parameter of EnforceRotCond is the index past actual end.
    int start = 0;
    int end = islower ? _grid_cond.dims[dir] - 1 : _grid_cond.dims[dir];
    EnforceRotCond( ft, Efield.data(dir), _E_cond.data(dir), start, end);
    EnforceRotCond( 1.0, Bfield.data(trans[0]), _B_cond.data(trans[0]), start, end);
    EnforceRotCond( 1.0, Bfield.data(trans[1]), _B_cond.data(trans[1]), start, end);

  }

  virtual FieldBCType bcType() const { return FieldBCType::ROTATING_CONDUCTOR; };

};
#endif // FIELDBC_ROTATING_CONDUCTOR_H
