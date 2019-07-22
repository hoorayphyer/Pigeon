#ifndef _IO_EXPORTEE_IMPL_HPP_
#define _IO_EXPORTEE_IMPL_HPP_

#include "io/exportee.hpp"

namespace io {
  void EB ( const mani::Grid<Real,DGrid>& grid; // local grid
            field::Field<Real, 3, DGrid>& E;
            field::Field<Real, 3, DGrid>& B;
            field::Field<RealJ, 3, DGrid>& J;// J is Jmesh on a replica
            particle::map<particle::array<Real,S>>& particles ) {
      field::Field<Real, 1, DGrid> tmp( E.mesh() );

      std::vector<int> dims(DGrid);
      for ( int i = 0; i < DGrid; ++i ) dims[i] = io_field.mesh().extent()[i];

      std::string varname;
      for ( int comp = 0; comp < 3; ++comp ) {
          downsample( DSRatio, io_field, E, comp );
          field::copy_sync_guard_cells( io_field, *cart_opt );
          varname = "E" + std::to_string(comp+1);
          pmpio([&](auto& dbfile){
                  dbfile.put_var( varname, MeshExport, io_field[0].data().data(), dims );
                });
          put_to_master( varname, cart_opt->size());
      }

      for ( int comp = 0; comp < 3; ++comp ) {
          downsample( DSRatio, io_field, B, comp );
          field::copy_sync_guard_cells( io_field, *cart_opt );
          varname = "B" + std::to_string(comp+1);
          pmpio([&](auto& dbfile){
                  dbfile.put_var( varname, MeshExport, io_field[0].data().data(), dims );
                });
          put_to_master( varname, cart_opt->size());
        }

        {
          for ( const auto& I : apt::Block(tmp.mesh().bulk_dims()) ) {
            tmp[0](I) = J[comp](I);
          }
          tmp.set_offset( 0, J[comp].offset() );
          { // normalize J to orthonormal basis
            // define a function pointer.
            Real(*h_func)(Real,Real,Real) = nullptr;
            switch(comp) {
            case 0: h_func = Metric::template hh<0,Real>; break;
            case 1: h_func = Metric::template hh<1,Real>; break;
            case 2: h_func = Metric::template hh<2,Real>; break;
            }

            static_assert( DGrid == 2 );
            const auto& ofs = tmp[0].offset();
            apt::array<Real,DGrid> q {};
            for ( const auto& I : apt::Block(tmp.mesh().bulk_dims()) ) {
              // TODOL use generator
              for ( int i = 0; i < DGrid; ++i ) q[i] = grid[i].absc(I[i], ofs[i]);

              auto h = h_func(q[0], q[1], 0.0);
              if ( std::abs(h) > 1e-12 )
                tmp[0](I) /= h;
              else
                tmp[0](I) = 0.0;
            }
          }
          downsample( DSRatio, io_field, tmp, 0 ); // 0 is to tmp
          field::copy_sync_guard_cells( io_field, *cart_opt );
          varname = "J"+std::to_string(comp+1);
          pmpio([&](auto& dbfile){
                  dbfile.put_var( varname, MeshExport, io_field[0].data().data(), dims );
                });
          put_to_master( varname, cart_opt->size());
        }
      }
  }
}

#endif
