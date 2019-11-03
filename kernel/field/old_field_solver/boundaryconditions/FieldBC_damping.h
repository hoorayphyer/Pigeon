#ifndef FIELDBC_DAMPING_H
#define FIELDBC_DAMPING_H

#include "FieldBC.h"
#include <cmath>

class FieldBC_damping : public FieldBC {
private:
  const BoundaryPosition _bpos;
  Scalar _absorb_rate;
  Grid _grid_damp;
  Index _shift; // store the index in the original grid that corresponds to new grid's (0,0,0)
  VectorField<Scalar> _B_bg; // initial background B
  VectorField<Scalar> _E_bg;

  inline void Damp( Scalar& f_target, const Scalar& f_bg, const Scalar lambda ) const {
    f_target = ( f_target - f_bg ) * lambda + f_bg;
  }

  inline Scalar damping_profile( Scalar x ) {
    return x * x / 2.0;
  }

  void InitBackground ( const VectorField<Scalar>& Efield_bg, const VectorField<Scalar>& Bfield_bg ) {
    for ( int k = 0; k < _grid_damp.dims[2]; ++k ) {
      for ( int j = 0; j < _grid_damp.dims[1]; ++j ) {
        for ( int i = 0; i < _grid_damp.dims[0]; ++i ) {
          int i_s = i + _shift[0];
          int j_s = j + _shift[1];
          int k_s = k + _shift[2];
          for ( int comp = 0; comp < VECTOR_DIM; ++comp ) {
            _E_bg( comp, i, j, k ) = Efield_bg( comp, i_s, j_s, k_s );
            _B_bg( comp, i, j, k ) = Bfield_bg( comp, i_s, j_s, k_s );
          }
        }
      }
    }
  }

  // TODO  TODO
  // void ReadFromSnapshot( const SaveSnapshotProxy& proxy ) {

  //   int size = _E_bg.grid().size();
  //   int size_snapshot = proxy.fBC_damping_E_bg.size() / 3;

  //   if ( size_snapshot != size ) {
  //     throw std::runtime_error( "array size mismatch when reading in damping background at boundary " + std::to_string(_bpos) + "! Local has size " + std::to_string(size) + ", whereas snapshot has size" + std::to_string(size_snapshot) );
  //   }

  //   // copy the data from proxy
  //   for ( unsigned int i = 0; i < size; ++i ) {
  //     for ( int j = 0; j < 3; ++j ) {
  //       _E_bg.data(j)[i] = proxy.fBC_damping_E_bg[ 3 * i + j ];
  //       _B_bg.data(j)[i] = proxy.fBC_damping_B_bg[ 3 * i + j ];
  //     }
  //   }

  //   Logger::print_screen( Logger::gVerbosityLvl, "Rank", Logger::thisRank, "reads in damping background at boundary", std::to_string(_bpos) );

  // }

public:
  template<typename DBPane>
  FieldBC_damping( BoundaryPosition bpos, const Grid& grid, const DBPane& pane, const VectorField<Scalar>& Efield_bg, const VectorField<Scalar>& Bfield_bg )
    : _bpos(bpos), _grid_damp(grid) {
    _absorb_rate = pane.damping_rate;

    int dir = _bpos / 2;
    bool isupper = _bpos % 2;

    // taylor the full grid to fit the damping region.
    // guard cells are irrelevant for _grid_damp. However, we still make
    // the grid_damp contain guard cells in the interested direction while
    // the grid dim is equal to indent.
    _grid_damp.dims[dir] = grid.indent[bpos];
    _grid_damp.indent[dir * 2] = _grid_damp.guard[dir];
    _grid_damp.indent[dir * 2 + 1] = _grid_damp.guard[dir];

    if( isupper ) {
        _grid_damp.lower[dir] = grid.pos( dir, grid.dims[dir] - grid.indent[bpos] - 1, 1 );
        _shift[dir] = grid.dims[dir] - grid.indent[bpos];
    }

    // resize bg fields and initialize
    _E_bg.resize( _grid_damp );
    _B_bg.resize( _grid_damp );

    // TODO TODO
    // const auto& proxy = InfoCollector::Instance().ssProxy;
    // if ( proxy.isReadDampingBg ) {
    //   ReadFromSnapshot(proxy);
    // } else {

    // TODO InitBackground is needed
      // InitBackground( Efield_bg, Bfield_bg );

    // }

  }

  ~FieldBC_damping() override = default;

  virtual void Apply(Scalar dt, VectorField<Scalar>& Efield, VectorField<Scalar>& Bfield, Scalar time) override {
    int dir = _bpos / 2;
    bool isupper = _bpos % 2;

    // TODO: the following implementation assumes upper boundary.
    Scalar r_out = std::exp( _grid_damp.pos(dir, _grid_damp.dims[dir] - 1, 1) );
    Scalar r_indent = std::exp( _grid_damp.lower[dir] );
    Scalar thickness = r_out - r_indent;
    int trans[2] = { (dir + 1) % VECTOR_DIM , (dir + 2) % VECTOR_DIM };

    // store the staggers of E and B along dir
    Index E_stags = { GetStagProperty( FieldType::ETYPE, 0 )[dir],
                      GetStagProperty( FieldType::ETYPE, 1 )[dir],
                      GetStagProperty( FieldType::ETYPE, 2 )[dir] };

    Index B_stags = { GetStagProperty( FieldType::BTYPE, 0 )[dir],
                      GetStagProperty( FieldType::BTYPE, 1 )[dir],
                      GetStagProperty( FieldType::BTYPE, 2 )[dir] };
    // lambda[0] corresponds to unstagger, lambda[1] to stagger.
    Scalar lambda[2] = { 0.0, 0.0 };

    for ( int i = 0; i < _grid_damp.dims[dir]; ++i ) {
      int i_s = i + _shift[dir];
      Scalar q1 = _grid_damp.pos(0, i, 0); // unstaggered q1
      Scalar q1_s = _grid_damp.pos(0, i, 1); // staggered q1
      Scalar h = std::exp( q1 ) - r_indent;
      Scalar h_s = std::exp( q1_s ) - r_indent;
      lambda[0] = 1.0 - _absorb_rate * dt * damping_profile( h / thickness );
      lambda[1] = 1.0 - _absorb_rate * dt * damping_profile( h_s / thickness );

      for ( int k = 0; k < _grid_damp.dims[trans[1]]; ++k ) {
        int k_s = k + _shift[trans[1]];
        for ( int j = 0; j < _grid_damp.dims[trans[0]]; ++j ) {
          int j_s = j + _shift[trans[0]];

          // damp transverse fields
          for ( int n = 0; n < VECTOR_DIM - 1; ++n ) {
            Damp( Efield( trans[n], i_s, j_s, k_s ), _E_bg( trans[n], i, j, k ),
                  lambda[ E_stags[ trans[n] ] ] );
            Damp( Bfield( trans[n], i_s, j_s, k_s ), _B_bg( trans[n], i, j, k ),
                  lambda[ B_stags[ trans[n] ] ] );
          }

          // damp normal fields
          Damp( Efield( dir, i_s, j_s, k_s ), _E_bg( dir, i, j, k ),
                lambda[ E_stags[dir] ] );
          Damp( Bfield( dir, i_s, j_s, k_s ), _B_bg( dir, i, j, k ),
                lambda[ B_stags[dir] ] );
        }
      }
    }

    // // don't push particles that already entered damping layer
    // const Grid& grid_full = data.Efield.grid();
    // for ( auto& elm : data.particles ) {
    //   auto& particles = elm.second;
    //   // ignore neutral particles
    //   if ( std::abs( particles.Attributes().charge ) < 1e-10 ) continue;

    //   for ( int idx = 0; idx < particles.Number(); ++idx ) {
    //     auto& ptc = particles.PtcData()[idx];
    //     const int& C_dir = grid_full.getCell( ptc.cell )[dir];
    //     // check if the particle lies inside the layer
    //     if ( C_dir >= _shift[dir] && C_dir < _shift[dir] + _grid_damp.dims[dir] ) {
    //       set_bit( ptc.flag, ParticleFlag::ignore_force );
    //     }
    //   }

    // }

  }

  // TODO TODO
  // For saving snapshots. This function is not a virtual
  // void CopyToSnapshot( SaveSnapshotProxy& proxy ) const {
  //   unsigned int size = _E_bg.grid().size();

  //   if ( size != proxy.fBC_damping_E_bg.size() ) {
  //     proxy.fBC_damping_E_bg.resize( 3 * size );
  //     proxy.fBC_damping_B_bg.resize( 3 * size );
  //   }

  //   // copy data into proxy
  //   for ( unsigned int i = 0; i < size; ++i ) {
  //     for( int j = 0; j < 3; ++j ) {
  //       proxy.fBC_damping_E_bg[ 3 * i + j ] = _E_bg.data(j)[i];
  //       proxy.fBC_damping_B_bg[ 3 * i + j ] = _B_bg.data(j)[i];
  //     }
  //   }

  // }

};


#endif // FIELDBC_DAMPING_H
