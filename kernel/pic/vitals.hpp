#ifndef _PIC_VITALS_HPP_
#define _PIC_VITALS_HPP_

#include "particle/properties.hpp"
#include "timer/timer.hpp"
// #include "msh/mesh_shape_interplay_impl.hpp" // WeightFinder // FIXME
#include <fstream>

// FIXME
// namespace pic {
//   template < int DGrid,
//              typename T,
//              template < typename > class S,
//              class ShapeF,
//              class Metric
//              >
//   std::optional<field::Field<T,1,DGrid>> skip_depth( const particle::map<particle::Properties>& properties,
//                                                      const particle::map<particle::array<T,S>>& particles,
//                                                      const dye::Ensemble<DGrid>& ens,
//                                                      const std::optional<mpi::CartComm>& cart,
//                                                      const apt::Grid<T,DGrid>& grid,
//                                                      T classic_electron_radius
//                   ) {
//     field::Field<T,1,DGrid> tmp;
//     apt::array< field::offset_t, DGrid > ofs;
//     apt::Index<DGrid> ext;
//     {
//       for ( int i = 0; i < DGrid; ++i ) ofs[i] = MIDWAY;
//       for ( int i = 0; i < DGrid; ++i ) ext[i] = ShapeF::support();
//       tmp.set_offset(0, ofs);
//     }

//     // calculate wp^2. Note there is 4\pi re in front
//     for ( auto sp : particles ) {
//       if ( properties[sp].charge_x == 0 ) continue;

//       T e2m = properties[sp].charge_x == 0;
//       e2m = e2m * e2m / properties[sp].mass_x;
//       e2m *= 4.0 * std::acos(-1.0l) * classic_electron_radius;

//       for ( const auto& ptc : particles[sp] ) {
//         msh::impl::WeightFinder<T, DGrid,S<T>::Dim, ShapeF> wf(msh::to_standard(grid, ptc.q()), ofs, ShapeF() );
//         for ( const auto& I : apt::Block(ext) ) {
//           apt::Index<DGrid> Ids;
//           for ( int i = 0; i < DGrid; ++i ) Ids[i] = wf.I_b()[i] + I[i];
//           tmp[0](Ids) += ptc.frac() * e2m * wf.weight(I);
//         }
//       }
//     }
//     // reduce to ens chief
//     ens.reduce_to_chief( mpi::by::SUM, tmp[0].data().data(), tmp[0].size() );
//     if ( !ens.is_chief() ) return {};
//     field::copy_sync_guard_cells( tmp, *cart );
//     // turn to skin depth. Include scale factors
//     auto dV = dV(grid);
//     for ( const auto& I : apt::Block(ext) ) {
//       apt::array<T,3> q {0,0,0};
//       for ( int i = 0; i < DGrid; ++i ) q[i] = grid[i].absc(I[i],ofs[i]);

//       // find the largest hdx
//       apt::array<T,3> hdx {0,0,0};
//       hdx[0] = Metric::h<0>( q[0], q[1], q[2] ) * grid[0].delta();
//       if constexpr ( DGrid > 1 ) {
//           hdx[1] = Metric::h<1>( q[0], q[1], q[2] ) * grid[1].delta();
//           if constexpr ( DGrid > 2 ) {
//               hdx[2] = Metric::h<2>( q[0], q[1], q[2] ) * grid[2].delta();
//             }
//         }
//       for ( int i = 1; i < 3; ++i  )
//         hdx[0] = std::max(hdx[0], hdx[i]);

//       tmp[0](I) /= (Metric::hhh(q[0], q[1], q[2]) * dV);

//     }
//     return {tmp};
//   }
// }

namespace pic {
  namespace vital {
    // significant only on world rank 0
    std::vector<double> num_ptcs_prev;
    std::vector<double> num_scat_prev;
    double t_phys_prev = 0;
  }

  template < int DGrid,
             typename T,
             template < typename > class S
             >
  void check_vitals( std::string filename, T t_phys, const dye::Ensemble<DGrid>& ens,
                     const std::optional<mpi::CartComm>& cart,
                     const particle::map<particle::Properties>& properties,
                     const particle::map<particle::array<T,S>>& particles,
                     const particle::map<T>& N_scat ) {
    using namespace particle;
    static int counter = 0;
    constexpr int interval = 40;
    static tmr::Timestamp stopwatch;

    std::vector<double> buffer;
    {
      for ( auto sp : particles ) {
        double num = 0.0;
        for ( const auto& ptc : particles[sp] ) {
          if ( ptc.is(flag::exist) ) num += ptc.frac();
        }
        buffer.push_back(num);
      }
      for ( auto sp : N_scat ) buffer.push_back(N_scat[sp]);
      ens.reduce_to_chief( mpi::by::SUM, buffer.data(), buffer.size() );
    }
    if ( !cart ) return;

    buffer.push_back(ens.size());
    cart->template reduce<true>( mpi::by::SUM, 0, buffer.data(), buffer.size() );

    if ( cart->rank() != 0 ) return;

    auto* p1 = buffer.data();
    auto* p2 = buffer.data() + particles.size();

    {
      std::ofstream out(filename, std::ios_base::app);
      if ( counter % interval == 0 ) {
        out << "t_phys|\tnprocs|\tlapse/hr|\tTotal load(rate)|\tcumulative new ptcs from scattering(rate)" << std::endl;
        out << "species ordering : ";
        for ( auto sp : N_scat )
          out << properties[sp].nickname << " ";
        out << std::endl;
        counter = 0;
      }
      ++counter;
      out << apt::fmt("%.2f",t_phys) << "|\t" << buffer.back() << "|\t" << apt::fmt("%8.2f", stopwatch.lapse().in_units_of("s").val() / 3600.0) <<  "|\t";
      for ( int i = 0; i < particles.size(); ++ i ) {
        out << apt::fmt("%.2e", p1[i]) << "(" << apt::fmt( "%.2e", (p1[i] - vital::num_ptcs_prev[i]) / (t_phys- vital::t_phys_prev) ) << ") ";
        vital::num_ptcs_prev[i] = p1[i];
      }
      out << "|\t";
      for ( int i = 0; i < N_scat.size(); ++i ) {
        out << apt::fmt("%.2e", p2[i]) << "(" << apt::fmt( "%.2e", (p2[i] - vital::num_scat_prev[i]) / (t_phys- vital::t_phys_prev) ) <<  ") ";
        vital::num_scat_prev[i] = p2[i];
      }
      out << std::endl;
      out.close();

      vital::t_phys_prev = t_phys;
    }

    return;
  }
}

#endif
