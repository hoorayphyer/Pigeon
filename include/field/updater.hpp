#ifndef  _FIELD_UPDATER_HPP_
#define  _FIELD_UPDATER_HPP_

#include "kernel/coordinate.hpp"


// TODO double check the new staggering
// TODO check 4\pi
namespace field {
  template < typename Real, knl::coordsys_t CS >
  class Updater {
  private:
    void ComputeBfieldUpdate( Scalar dt, VectorField<Scalar>& Efield, VectorField<Scalar>& Bfield, const VectorField<Scalar>& current) {
      // FIXME: There is something wrong with the extent of the
      // laplacian. Maybe I neglected the distinction between staggered
      // and unstaggered fields

      _BfieldOld.copyFrom(Bfield);

      _domain->SendGuardCells(_BfieldOld);

      if (std::abs(_beta) > EPS) {
        _fd.ComputeLaplacian(Bfield, _BfieldDelta, FieldType::BTYPE, _is_at_boundary.data(), *_domain, _d_start, _d_ext, true);

        _BfieldDelta.multiplyBy(_alpha * _beta * dt * dt);

        Bfield.addBy(_BfieldDelta);
      }

      _domain->SendGuardCells(Efield);

      // FIXME: using _isBdry_EJ as an expediency. Should use _domain->isBoundary().
      bool isBdry_E[NUM_BOUNDARIES] = {};
      for ( int b = 0; b < NUM_BOUNDARIES; ++b ) {
        isBdry_E[b] = _is_at_boundary[b];
      }

      if ( _is_at_boundary[0] ) {
        isBdry_E[0] = false;
        _d_start[0] -= 1;
        _d_ext[0] += 1;
      }
      _fd.ComputeCurl(Efield, _BfieldDelta, FieldType::ETYPE, isBdry_E, _d_start, _d_ext);
      _BfieldDelta.multiplyBy(-1.0 * dt);
      Bfield.addBy(_BfieldDelta);


      if ( _is_at_boundary[0] ) {
        _d_start[0] += 1;
        _d_ext[0] -= 1;
        isBdry_E[0] = true;
      }

      _fd.ComputeCurl(current, _BfieldDelta, FieldType::ETYPE, isBdry_E, _d_start, _d_ext);
      _BfieldDelta.multiplyBy(_alpha * dt * dt);
      Bfield.addBy(_BfieldDelta);

      _domain->SendGuardCells(Bfield);

      _BfieldDelta.copyFrom(Bfield);

      Scalar factor = dt * dt * _alpha * _alpha;
      for (int j = 0; j < 5; ++j) {
        _BfieldDelta.multiplyBy(factor);
        _fd.ComputeLaplacian(_BfieldDelta, _BfieldDelta, FieldType::BTYPE, _is_at_boundary.data(), *_domain, _d_start, _d_ext, true);

        Bfield.addBy(_BfieldDelta);
        _domain->SendGuardCells(_BfieldDelta);
      }

      _domain->SendGuardCells(Bfield);
    }

    void ComputeEfieldUpdate(Scalar dt, VectorField<Scalar>& Bfield, VectorField<Scalar>& Efield, const VectorField<Scalar>& current) {
      _EfieldDelta.assign(0.0);
      //  _domain->SendGuardCells(Bfield);

      _fd.ComputeCurl(Bfield, _EfieldTmp, FieldType::BTYPE, _is_at_boundary.data(), _d_start, _d_ext);
      _EfieldTmp.multiplyBy(_alpha);
      _EfieldDelta.addBy(_EfieldTmp);

      _fd.ComputeCurl(_BfieldOld, _EfieldTmp, FieldType::BTYPE, _is_at_boundary.data(), _d_start, _d_ext);
      _EfieldTmp.multiplyBy(_beta);
      _EfieldDelta.addBy(_EfieldTmp);

      _EfieldDelta.subtractBy(current);
      _EfieldDelta.multiplyBy(dt);

      Efield.addBy(_EfieldDelta);

      for (int i = 0; i < VECTOR_DIM; i++)
        CleanField(Efield.data(i), UPPER_LIMIT);

      _domain->SendGuardCells(Efield);
    }

    void SetFieldBC(std::unordered_map<BoundaryPosition, FieldBC*>& field_bdry_cond, const AperParams& params, Domain* domain, const DBPane_BoundaryCondition& pane ) {
      const auto& isBoundary = params.ens_specs.is_at_boundary;
      for ( const auto& elm : pane.fieldBC ) {
        auto b = elm.first;
        if ( !isBoundary[b] ) continue;

        const auto& fbc = elm.second;
        switch (fbc.type) {
        case FieldBCType::COORDINATE :
          field_bdry_cond[b] = new FieldBC_coordinate(b);
          break;
        case FieldBCType::DAMPING :
          field_bdry_cond[b] = new FieldBC_damping( b, params.dt, params.grid, pane.fieldBC.at(b), VectorField<Scalar>(), VectorField<Scalar>() );
          break;
        case FieldBCType::ROTATING_CONDUCTOR :
          field_bdry_cond[b] = new FieldBC_rotating_conductor( b, params, fbc );
          break;
        case FieldBCType::DIFFERENTIAL_ROTATION :
          // FIXME fix interface?
          // field_bdry_cond[b] = new FieldBC_differential_rotation( b, params, db.coordType,
          //                                                  fbc.bcFunc.ft, fbc.bcFunc.B1,
          //                                                  fbc.bcFunc.B2, fbc.bcFunc.B3 );
          break;
        case FieldBCType::FUNCTION_SPLIT :
          // FIXME what is this?
          // field_bdry_cond = new FieldBC_functionsplit( b, params );
          break;
        case FieldBCType::PML :
          // FIXME update interface
          // field_bdry_cond = new FieldBC_pml( b, domain, params, scales );
          break;
        } // end of switch
      } // end of for
    }

  public:
    void operator( auto& E, auto& B, const auto& J, Real dt ) () {

      ComputeBfieldUpdate(dt, E, B, J);
      ComputeEfieldUpdate(dt, B, E, J);
      // TODO
      // for ( auto& elm : _fieldBC )
      //   elm.second->Apply( data, timestep * dt );
    }

  };
}

#endif
