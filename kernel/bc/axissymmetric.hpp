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
      // TODO Guard cells values are needed when doing interpolating E and B
      // E_theta, B_r, B_phi are on the axis. All but B_r should be set to zero
      const auto& mesh = _Efield.mesh();
      if ( _is_at_axis_lower ) {
        int n = 0;
        for ( const auto& trI : mesh.project(1, {}, mesh.extent() ) ) {
          _Efield[1][trI | n] = 0.0;
          _Bfield[2][trI | n] = 0.0;
        }
      }

      if ( _is_at_axis_upper ) {
        int n = _grid[1].dim();
        for ( const auto& trI : mesh.project(1, {}, mesh.extent() ) ) {
          _Efield[1][trI | n] = 0.0;
          _Bfield[2][trI | n] = 0.0;
        }
      }
    }
  };
}

#endif
