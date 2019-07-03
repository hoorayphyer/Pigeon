#ifndef _BC_AXISSYMMETRIC_HPP_
#define _BC_AXISSYMMETRIC_HPP_

namespace bc {
  template < int DGrid,
             typename Real,
             template < typename > class Specs,
             typename RealJ >
  struct Axissymmetric {
  private:
    const mani::Grid<Real,DGrid>& _grid;
    field::Field<Real, 3, DGrid>& _Efield;
    field::Field<Real, 3, DGrid>& _Bfield;

    const int axis_dir = 1;
    bool _is_at_axis_lower = false;
    bool _is_at_axis_upper = false;

  public:
    Axissymmetric ( const mani::Grid<Real,DGrid>& localgrid,
                    field::Field<Real, 3, DGrid>& Efield,
                    field::Field<Real, 3, DGrid>& Bfield,
                    const field::Field<RealJ, 3, DGrid>& Jfield,
                    const particle::map<particle::array<Real,Specs>>& particles )
      : _grid(localgrid), _Efield(Efield), _Bfield(Bfield) {
      _is_at_axis_lower = std::abs( localgrid[axis_dir].lower() - 0.0 ) < localgrid[axis_dir].delta();
      _is_at_axis_upper = std::abs( localgrid[axis_dir].upper() - std::acos(-1.0) ) < localgrid[axis_dir].delta();
    }

    void operator() () {
      // NOTE Guard cells values are needed when interpolating E and B
      // E_theta, B_r, B_phi are on the axis. All but B_r should be set to zero
      const auto& mesh = _Efield.mesh();
      if ( _is_at_axis_lower ) {
        for ( const auto& trI : mesh.project(axis_dir, mesh.origin(), mesh.extent() ) ) {
          for ( int n = mesh.origin()[axis_dir]; n < 0; ++n ) {
            // MIDWAY in axis_dir
            _Efield[0][trI | n] = _Efield[0][ trI | ( -n - 1 ) ];
            _Efield[2][trI | n] = - _Efield[2][ trI | ( -n - 1 ) ];
            _Bfield[1][trI | n] = - _Bfield[1][ trI | ( -n - 1 ) ];

            // INSITU in axis_dir
            _Efield[1][trI | n] = - _Efield[1][ trI | -n ];
            _Bfield[0][trI | n] = _Bfield[0][ trI | -n ];
            _Bfield[2][trI | n] = - _Bfield[2][ trI | -n ];
          }

          _Efield[1][trI | 0] = 0.0;
          _Bfield[2][trI | 0] = 0.0;
        }
      }

      if ( _is_at_axis_upper ) {
        int bulk = mesh.bulk_dim(axis_dir);
        for ( const auto& trI : mesh.project(axis_dir, mesh.origin(), mesh.extent() ) ) {
          for ( int n = bulk; n < bulk + mesh.guard(); ++n ) {
            // MIDWAY in axis_dir
            _Efield[0][trI | n] = _Efield[0][ trI | ( 2*bulk - 1 - n ) ];
            _Efield[2][trI | n] = - _Efield[2][ trI | ( 2*bulk - 1 - n ) ];
            _Bfield[1][trI | n] = - _Bfield[1][ trI | ( 2*bulk - 1 - n ) ];

            // INSITU in axis_dir
            _Efield[1][trI | n] = - _Efield[1][ trI | ( 2*bulk -n ) ];
            _Bfield[0][trI | n] = _Bfield[0][ trI | ( 2*bulk -n ) ];
            _Bfield[2][trI | n] = - _Bfield[2][ trI | ( 2*bulk -n ) ];
          }

          _Efield[1][trI | bulk] = 0.0;
          _Bfield[2][trI | bulk] = 0.0;
        }
      }
    }
  };
}

#endif
