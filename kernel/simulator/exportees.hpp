#ifndef _APERTURE_EXPORTEES_HPP_
#define _APERTURE_EXPORTEES_HPP_

#include "io/exportee.hpp"
#include "particle/array.hpp"

namespace aperture {
  struct FldComp : FieldBasedExportee {
  private:

  public:
    FldComp( const field::Field& fld, int comp = 0 );
    // TODO do downsampling here
    T operator() ( ) const override {
      
    }
  };

  // TODO average to expf should factor in the scale functions, i.e. one should find the downsampled value by conserving the flux.
  // fd.ComputeDivergence( Efield, _scalar_tmp, FieldType::ETYPE, params.ens_specs.is_bdry.data() );
  // domain.SendGuardCells( _scalar_tmp );
  // put_var(dbfile, "divE", _scalar_tmp.ptr());

  // fd.ComputeDivergence( Bfield, _scalar_tmp, FieldType::BTYPE, params.ens_specs.is_bdry.data() );
  // domain.SendGuardCells( _scalar_tmp );
  // put_var(dbfile, "divB", _scalar_tmp.ptr());

  // GetFlux
}

namespace particle {
  template < typename T, int DPtc, typename state_t, int DGrid >
  struct ArrayExportee : io::ParticleBasedExportee<T,DGrid> {
  protected:
    const array<T,DPtc,state_t>& _array;

  public:
    constexpr ParticleArrayExportee( const array<T,DPtc,state_t>& array )
      : _array(array) {}

    load_t begin() const override final {
      load_t res = 0;
      return _array[res].is(flag::empty) ? next(res) : res;
    }

    load_t end() const override final { return _array.size(); }

    load_t& next( load_t& index ) const {
      ++index;
      while ( _array[index].is(flag::empty) && index != end() ) ++index;
      return index;
    }

    apt::array<T, DGrid> loc() ( load_t index ) const override final {
      apt::array<T, DGrid> res;
      apt::foreach<0,DGrid>
        ( []( auto& r, auto q ) { r = q;}, res, _array[index].q() );
      return res;
    }

  };
}

namespace aperture {

  // TODO check if the species is actually used. NOTE this is not equivalent to having zero particles
  // TODO RescaledShapeF to be associated with DataExporter Machine, or export_mesh
  template < typename ShapeF >
  struct RescaledShapeF {
    RescaledShapeF shapef() const {};
  };

  template < typename T >
  using ArrayExportee = particle::ArrayExportee<T, traits::DPtc, traits::state_t, traits::DGrid >;


  template < typename T >
  struct Number : particle::ArrayExportee<T> {
    using particle::ArrayExportee<T>::ArrayExportee<T>;
    T calc() const override {  return 1.0; }
  };

  template < typename T >
  struct Energy : particle::ArrayExportee<T> {
  public:
    using particle::ArrayExportee<T>::ArrayExportee<T>;
    // EnergyDensity( particle::species sp );
    T calc( const auto& ptc ) const override {
      return std::sqrt( ( mass(sp) > 0 ) + apt::sqabs(ptc.p()) );
    }
  };

  template < typename T >
  struct Momentum : particle::ArrayExportee<T> {
  public:
    using particle::ArrayExportee<T>::ArrayExportee<T>;
    Momentum( int comp );

    T calc( const auto& ptc ) const override {
      return ptc.p()[_comp];
    }
  };

  // TODO shapef should be with respect to datagrid. Rescaling?
  // FIXME TODO save gamma P density, which is total gamma / physical cell volume. Later this divided by number density gives us gamma P per unit particle. TODO: think this over.

  // TODO need send add cells
  // TODO boundary conditions?? YES
  // FIXME TODO save number density instead of just number

        // // charged species do the following. Check charge
        // // FIXME in the presence of this, may we just drop exporting J_total?
        // if ( is_charged(ptcType) ) {
        //   auto& J_sp = _vector_tmp;
        //   J_sp.assign(0.0);
        //   CurrentDepositer j_dep;
        //   std::visit([&]( auto& x ) {
        //                j_dep( J_sp, dt, x, params, comm, domain );
        //              }, ptcs );
        //   if ( _comm -> is_primary() ) {
        //     put_var( dbfile, "J1_"+ptcStr, J_sp.ptr(0) );
        //     put_var( dbfile, "J2_"+ptcStr, J_sp.ptr(1) );
        //     put_var( dbfile, "J3_"+ptcStr, J_sp.ptr(2) );
        //   }

        // }

        // { auto& pc_rate = _scalar_tmp;
        //   pc_rate.copyFrom(data.pairCreationEvents);
        //   domain.EnsReduceFields(pc_rate);
        //   if ( _comm -> is_primary() ) {
        //     // convert to rate
        //     pc_rate.multiplyBy( 1.0 / ( _pane.interval * params.dt ) );
        //     domain.SendGuardCells( pc_rate );
        //     put_var( dbfile, "Pair_Creation_Rate", pc_rate.ptr() );
        //   }
        // }
}

#endif
