#ifndef _BC_AXISSYMMETRIC_HPP_
#define _BC_AXISSYMMETRIC_HPP_

namespace bc {
  template < int DGrid,
             typename Real,
             template < typename > class Specs,
             typename RealJ >
  struct Axissymmetric {
  private:
    static_assert(DGrid==2);
    const mani::Grid<Real,DGrid>& _grid;
    field::Field<Real, 3, DGrid>& _E;
    field::Field<Real, 3, DGrid>& _B;
    field::Field<RealJ, 3, DGrid>& _J;

    static constexpr int AxisDir = 1;
    bool _is_lower = false;
    bool _is_upper = false;

  public:
    // NOTE f is a function pointer to void (*)( T&, T& ). The complicated expression is just to make template deduction work.
    template < typename T >
    static void symmetrize( bool is_at_axis_lower, bool is_at_axis_upper,
                            field::Component<T,DGrid,false> comp, // TODOL semantics on comp
                            void (*f)( decltype(comp[0])& val_guard, decltype(comp[0])& val_bulk ) ) {
      const auto& mesh = comp.mesh();
      int mirror_sum = (comp.offset()[AxisDir] == MIDWAY ) ? -1 : 0;

      if ( is_at_axis_lower ) {
        for ( const auto& trI : mesh.project(AxisDir, mesh.origin(), mesh.extent() ) ) {
          for ( int n = mesh.origin()[AxisDir]; n < mirror_sum + 1; ++n ) { // n < 0 if MIDWAY, n < 1 is INSITU
            f( comp[ trI | n ], comp[ trI | (mirror_sum - n) ] );
          }
        }
      }
      if ( is_at_axis_upper ) {
        const int bulk = mesh.bulk_dim(AxisDir);
        for ( const auto& trI : mesh.project(AxisDir, mesh.origin(), mesh.extent() ) ) {
          for ( int n = bulk; n < bulk + mesh.guard(); ++n ) {
            f( comp[ trI | n ], comp[ trI | 2*bulk + mirror_sum - n ] );
          }
        }
      }
    }


    Axissymmetric ( const mani::Grid<Real,DGrid>& localgrid,
                    field::Field<Real, 3, DGrid>& Efield,
                    field::Field<Real, 3, DGrid>& Bfield,
                    field::Field<RealJ, 3, DGrid>& Jfield,
                    const particle::map<particle::array<Real,Specs>>& particles )
      : _grid(localgrid), _E(Efield), _B(Bfield), _J(Jfield) {
      _is_lower = std::abs( localgrid[AxisDir].lower() - 0.0 ) < localgrid[AxisDir].delta();
      _is_upper = std::abs( localgrid[AxisDir].upper() - std::acos(-1.0) ) < localgrid[AxisDir].delta();
    }

    void setEB () {
      // NOTE Guard cells values are needed when interpolating E and B
      // E_theta, B_r, B_phi are on the axis. All but B_r should be set to zero
      auto assign = []( Real& v_g, Real& v_b ) noexcept { v_g = v_b; };
      auto neg_assign = []( Real& v_g, Real& v_b ) noexcept {
                          v_g = ( &v_g == &v_b ) ? 0.0 : - v_b;
                        };
      // MIDWAY in AxisDir
      symmetrize(_is_lower, _is_upper, _E[0], assign );
      symmetrize(_is_lower, _is_upper, _E[2], neg_assign );
      symmetrize(_is_lower, _is_upper, _B[1], neg_assign );

      // INSITU in AxisDir
      symmetrize(_is_lower, _is_upper, _E[1], neg_assign );
      symmetrize(_is_lower, _is_upper, _B[0], assign );
      symmetrize(_is_lower, _is_upper, _B[2], neg_assign );
    }

    void setJ() {
      auto add_assign =
        []( RealJ& a, RealJ& b ) noexcept {
          a += b;
          b = a;
        };

      auto sub_assign =
        []( RealJ& a, RealJ& b ) noexcept {
          a -= b;
          b = -a;
        };
      // MIDWAY in AxisDir
      symmetrize(_is_lower, _is_upper, _J[0], add_assign );
      symmetrize(_is_lower, _is_upper, _J[2], sub_assign );
      // INSITU in AxisDir
      symmetrize(_is_lower, _is_upper, _J[1], sub_assign );
    }
  };
}

#endif
