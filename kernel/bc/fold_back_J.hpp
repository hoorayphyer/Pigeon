#ifndef _BC_FOLDBACKJ_HPP_
#define _BC_FOLDBACKJ_HPP_

namespace bc {
  template < int DGrid,
             typename Real,
             template < typename > class Specs,
             typename RealJ >
  struct FoldBackJ {
  private:
    static_assert(DGrid==2);
    const mani::Grid<Real,DGrid>& _grid;
    field::Field<RealJ, 3, DGrid>& _Jmesh;

    const int axis_dir = 1;
    bool _is_at_axis_lower = false;
    bool _is_at_axis_upper = false;

  public:
    FoldBackJ ( const mani::Grid<Real,DGrid>& localgrid,
                const field::Field<Real, 3, DGrid>& Efield,
                const field::Field<Real, 3, DGrid>& Bfield,
                field::Field<RealJ, 3, DGrid>& Jfield,
                const particle::map<particle::array<Real,Specs>>& particles )
      : _grid(localgrid), _Jmesh(Jfield) {
      _is_at_axis_lower = std::abs( localgrid[axis_dir].lower() - 0.0 ) < localgrid[axis_dir].delta();
      _is_at_axis_upper = std::abs( localgrid[axis_dir].upper() - std::acos(-1.0) ) < localgrid[axis_dir].delta();
    }

    void operator() () {
      const auto& mesh = _Jmesh.mesh();
      const int guard = mesh.guard();

      if ( _is_at_axis_lower ) {
        for ( const auto& trI : mesh.project(axis_dir, mesh.origin(), mesh.extent() ) ) {
          for ( int n = 0; n < guard; ++n ) {
            _Jmesh[0][ trI | n ] += _Jmesh[0][ trI | -1 - n ];
            _Jmesh[2][ trI | n ] -= _Jmesh[0][ trI | -1 - n ];
          }

          _Jmesh[1][trI | 0] = 0.0;
          // TODO check the negative sign here
          for ( int n = 0; n < guard; ++n )
            _Jmesh[1][ trI | 1 + n ] -= _Jmesh[1][ trI | -1 - n ];
        }
      }

      if ( _is_at_axis_upper ) {
        const int dim = mesh.bulk_dim(1);
        for ( const auto& trI : mesh.project(axis_dir, mesh.origin(), mesh.extent() ) ) {
          for ( int n = 0; n < guard; ++n ) {
            _Jmesh[0][ trI | dim - 1 - n ] += _Jmesh[0][ trI | dim + n ];
            _Jmesh[2][ trI | dim - 1 - n ] -= _Jmesh[0][ trI | dim + n ];
          }

          _Jmesh[1][trI | dim] = 0.0;
          for ( int n = 0; n < guard - 1; ++n ) // NOTE -1 is because of the favoring-lower convention
            _Jmesh[1][ trI | dim-1-n ] -= _Jmesh[1][ trI | dim+1+n ];

        }
      }
    }
  };
}

#endif
