#include "testfw/testfw.hpp"
#include "particle/annihilation_impl.hpp"
#include "particle/shapef.hpp"
#include "field/yee.hpp"
#include "msh/current_deposition_impl.hpp"

using namespace particle;
using R = float;
using aio::Specs;
using Vec_t = apt::Vec<R,Specs<R>::Dim>;

SCENARIO("Test Annihilation" , "[particle][mpi]") {
  const int num_procs = mpi::world.size();
  auto rwld_opt = aio::reduced_world(num_procs, mpi::world);
  if ( rwld_opt ) {
    constexpr int DGrid = 2;
    using RJ = R;

    const std::array<int,DGrid> dims {8,8};

    apt::Grid<R,DGrid> grid {{{0.0,1.0,dims[0]}, {0.0,1.0,dims[0]}}};

    field::Mesh<DGrid> mesh{{{ {0,dims[0],1},{0,dims[1],1} }}};
    field::Field<RJ,3,DGrid> J(mesh);
    for( int i = 0; i < 3; ++i ) {
      J.set_offset( i, field::yee::ofs_gen<DGrid>( field::yee::Etype, i ) );
    }

    auto init_ptcs =
      [&grid]( auto& ptcs, const auto frac, const auto Npc, const auto rel ) {
        for ( int j = 0; j < grid[1].dim(); ++j ) {
          R y = grid[1].absc(j,rel[1]);
          for ( int i = 0; i < grid[0].dim(); ++i ) {
            R x = grid[0].absc(i,rel[0]);
            for ( int c = 0; c < Npc; ++c ) {
              ptcs.emplace_back({x,y,0},{},frac);
            }
          }
        }
      };

    auto get_Npc_el
      = []( bool is_odd_rank ) {
          return is_odd_rank ? 2 : 10;
        };
    auto get_Npc_po
      = []( bool is_odd_rank ) {
          return is_odd_rank ? 6 : 1;
        };

    bool is_odd_rank = (rwld_opt->rank() % 2 == 1);

    array<R,Specs> el;
    R frac_el = 0.93;
    int Npc_el = get_Npc_el(is_odd_rank);
    std::array<R,DGrid> rel_el { 0.37, 0.53 };
    init_ptcs(el, frac_el, Npc_el, rel_el);

    array<R,Specs> po;
    R frac_po = 1.47;
    int Npc_po = get_Npc_po(is_odd_rank);
    std::array<R,DGrid> rel_po { 0.71, 0.19 };
    init_ptcs(po, frac_po, Npc_po, rel_po);

    R dt = 1.0;

    using ShapeF = shapef_t<shape::Cloud_In_Cell>;

    auto policy = []( auto n_e, auto n_p ) { return std::min(n_e,n_p); };

    annihilate(el,po,J,R(-1.0),R(1.0),grid,*rwld_opt,dt,ShapeF(), policy);

    auto count =
      [&mesh,&grid]( const auto& ptcs, bool ignore_frac = false ) {
        auto get_idx =
          [&] ( const auto& ptc ) {
            apt::Index<DGrid> I;
            for ( int i = 0; i < DGrid; ++i )
              I[i] = ( ptc.q(i) - grid[i].lower() ) / grid[i].delta();
            return mesh.linear_index(I);
          };
        field::Field<R,1,DGrid> num(mesh);
        if ( ignore_frac ) {
          for ( const auto& x : ptcs) {
            if ( not x.is(flag::exist) ) continue;
            num[0][get_idx(x)] += 1.0;
          }
        } else {
          for ( const auto& x : ptcs) {
            if ( not x.is(flag::exist) ) continue;
            num[0][get_idx(x)] += x.frac();
          }
        }
        return num;
      };

    { // test
      auto num_e = count(el);
      auto num_p = count(po);
      rwld_opt->reduce<true>(mpi::by::SUM, 0, num_e[0].data().data(), num_e[0].data().size() );
      rwld_opt->reduce<true>(mpi::by::SUM, 0, num_p[0].data().data(), num_p[0].data().size() );

      for ( int c = 0; c < 3; ++c )
        rwld_opt->reduce<true>(mpi::by::SUM, 0, J[c].data().data(), J[c].data().size() );

      R Npc_tot_el = frac_el * Npc_el;
      R Npc_tot_po = frac_po * Npc_po;

      rwld_opt->reduce<true>(mpi::by::SUM, 0, &Npc_tot_el, 1 );
      rwld_opt->reduce<true>(mpi::by::SUM, 0, &Npc_tot_po, 1 );

      auto num_e_nofr = count(el,true);
      auto num_p_nofr = count(po,true);
      rwld_opt->reduce<true>(mpi::by::SUM, 0, num_e_nofr[0].data().data(), num_e_nofr[0].data().size() );
      rwld_opt->reduce<true>(mpi::by::SUM, 0, num_p_nofr[0].data().data(), num_p_nofr[0].data().size() );

      if (rwld_opt->rank() == 0) {
        int size = rwld_opt->size();
        int num_odd = size / 2;
        int num_even = size - num_odd;
        REQUIRE( Npc_tot_el ==
                 Approx(frac_el *( get_Npc_el(true) * num_odd + get_Npc_el(false) * num_even ) ) );
        REQUIRE( Npc_tot_po ==
                 Approx(frac_po *( get_Npc_po(true) * num_odd + get_Npc_po(false) * num_even ) ) );

        R Nannih_pc = std::min(Npc_tot_el, Npc_tot_po);

        for ( int j = 0; j < grid[1].dim(); ++j ) {
          for ( int i = 0; i < grid[0].dim(); ++i ) {
            CHECK(num_e[0]({i,j}) == Approx(Npc_tot_el - Nannih_pc));
            CHECK(num_p[0]({i,j}) == Approx(Npc_tot_po - Nannih_pc));
          }
        }

        { // test count of objects, which matters more for performance
          R n_e {}, n_p {};
          if ( Npc_tot_el > Npc_tot_po ) {
            n_p = 0;
            n_e = std::ceil((Npc_tot_el - Nannih_pc) / frac_el);
          } else {
            n_e = 0;
            n_p = std::ceil((Npc_tot_po - Nannih_pc) / frac_po);
          }

          for ( int j = 0; j < grid[1].dim(); ++j ) {
            for ( int i = 0; i < grid[0].dim(); ++i ) {
              CHECK(num_e_nofr[0]({i,j}) == Approx(n_e));
              CHECK(num_p_nofr[0]({i,j}) == Approx(n_p));
            }
          }
        }

        { // test J
          Vec_t J_bulk;
          {
            auto Jtmp = J;
            Jtmp.reset();
            int I[2] = {mesh.range(0).full_size() / 2, mesh.range(1).full_size() / 2}; // get the center cell. Use full_size because q_std requires so
            apt::array<R,3> q1_std { I[0] + 0.5, I[1] + 0.5, 0.5 };
            apt::array<R,3> q0_std { I[0] +rel_el[0], I[1] + rel_el[1], q1_std[2] };
            msh::deposit( Jtmp, - R(1.0) / dt, ShapeF(), q0_std, q1_std );

            q0_std = { I[0] +rel_po[0], I[1] + rel_po[1], q1_std[2] };
            msh::deposit( Jtmp, R(1.0) / dt, ShapeF(), q0_std, q1_std );
            for ( int c = 0; c < 3; ++c )
              J_bulk[c] = std::accumulate(Jtmp[c].data().begin(),Jtmp[c].data().end(), R(0) );
            J_bulk *= Nannih_pc;
          }

          for ( int c = 0; c < 3; ++c ) {
            for ( int j = mesh.range(1).margin(0); j < mesh.range(1).size() - mesh.range(1).margin(1); ++j ) {
              for ( int i = mesh.range(0).margin(0); i < mesh.range(0).size() - mesh.range(0).margin(1); ++i ) {
                CAPTURE(c,i,j);
                CHECK( J[c]({i,j}) == Approx(J_bulk[c]) );
              }
            }
          }

        }
      }

    }

  }
  mpi::world.barrier();
}
