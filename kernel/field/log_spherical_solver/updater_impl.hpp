#include "field/log_spherical_solver/updater.hpp"
#include "field/sync.hpp" // FIXME not in use right now
#include "field/yee.hpp"
#include <optional>
#include <cassert>
#include <tuple>
#include <iostream> // FIXME

namespace field {
  constexpr int r_ = 0, th_ = 1, phi_ = 2;
  constexpr int bdry = 1; // set by 2nd order scheme

  template <typename R, int D>
  class _helper {
  private:
    mutable std::vector<R> _fr{};
    std::vector<R> _sin_th{};
    const int _b_r;
    const int _b_th;
    const R _estag; // exp( - 0.5 * dlnr )

  public:
    _helper(const apt::Grid<R, D> &grid, const int b_r, const int e_r,
            const int b_th, const int e_th)
        : _b_r(b_r), _b_th(b_th), _estag(std::exp(-grid[r_].delta()*0.5)) {
      _fr.resize(e_r - b_r);
      _fr.shrink_to_fit();
      _sin_th.resize(2 * (e_th - b_th));
      _sin_th.shrink_to_fit();

      for (int j = 0; j < e_th - b_th; ++j) {
        _sin_th[2 * j] = std::sin(grid[th_].absc(b_th + j));
        _sin_th[2 * j + 1] = std::sin(grid[th_].absc(b_th + j, 0.5));
      }
    }

    std::vector<R> &fr() const { return _fr; } // used mutable _fr
    auto &fr(int i) const { return _fr[i - _b_r]; } // used mutable _fr

    const auto &sin(int j, bool midway = false) const {
      return _sin_th[((j - _b_th) << 1) + midway];
    }
    const auto& estag() const {return _estag;}
  };


  template < typename R, int D >
  std::optional<const _helper<R,D>> hopt;

  template < typename R, int D >
  struct Diff_in_CurlE {
    const int mode; // 0 : d_r 1 : d_th
    const std::optional<int>& lo;
    const std::optional<int>& hi;
    const R dlnr;
    const R dth;

    using proxy_t =
        std::tuple<const int, const std::optional<int>&, const std::optional<int>&,
                   const Field<R,3,D> &, const int, R,
                   const int (&)[D], const int (&)[D], const bool >;
    static constexpr int Prefactor = 5;

    proxy_t operator()(const Field<R, 3, D> &E, const int comp, R prefactor,
                       const int (&b)[D], const int (&e)[D],
                       const bool continuous = false) {
      assert(comp != mode);
      prefactor /= ( mode == r_ ? dlnr : dth );
      return std::forward_as_tuple(mode,lo,hi,E,comp, prefactor,b,e,continuous);
    }

    static void d_r_dE(Component<R,D,false>&& output, proxy_t&& diff) {// TODOL sematics
      const auto &&[mode, lo, hi, E, comp, prefactor, b_in, e_in, continuous] = std::move(diff);
      const auto& m = E.mesh();
      int b = b_in[r_];
      int e = e_in[r_];

      if(!lo) b++;
      else if (not continuous) {
        R f = prefactor;
        b = *lo + bdry;

        for ( int j = b_in[th_]; j < e_in[th_]; ++j ) {
          int li_j = m.linear_index(th_, j);
          for ( int i = *lo; i < b; ++i ) {
            int li = li_j + (i+1) * m.stride(r_);
            // The following is O(dr^3) precise.
            auto x01 = E[comp][li - m.stride(r_)] - E[comp][li];
            auto x12 = E[comp][li] - E[comp][li + m.stride(r_)];
            auto x23 = E[comp][li + m.stride(r_)] - E[comp][li + (m.stride(r_) << 1)];
            output[li] += f * ( static_cast<R>(-71.0/24)*x01 + static_cast<R>(70.0/24)*x12 - static_cast<R>(23.0/24)*x23  );
          }
        }
      }

      {
        R f = prefactor;
        for (int j = b_in[th_]; j < e_in[th_]; ++j) {
          int li_j = m.linear_index(th_, j);
          for (int i = b; i < e; ++i) {
            int li = li_j + i * m.stride(r_);
            output[li] += f * (E[comp][li] - E[comp][li - m.stride(r_)]);
          }
        }
      }

      // NOTE at outer boundary, the r derivative on E_\theta and E_\phi is taken to vanish.
    }

    static void d_th_dE(Component<R,D,false>&& output, proxy_t&& diff) { // TODOL sematics
      const auto &&[mode, lo, hi, E, comp, prefactor, b_in, e_in, continuous] = std::move(diff);

      const auto& m = E.mesh();
      int b = b_in[th_];
      int e = e_in[th_];

      if(!lo) b++;
      else {
        b = bdry; // NOTE true for both comp == phi_ and comp == r_
        if (comp == phi_) {
          R f = prefactor / 3.0;
          for (int j = 0; j < b; ++j) {
            int li_j = m.linear_index(th_, j);
            for (int i = b_in[r_]; i < e_in[r_]; ++i) {
              int li = li_j + i * m.stride(r_);
              auto x = E[comp][li] + E[comp][li];
              x += x;
              x += x;
              x += E[comp][li];
              output[li] += f*( x - E[comp][li + m.stride(th_)]);
            }
          }
        } // NOTE if comp == r_, by symmetry output += 0, which doesn't need any action.
      }

      {
        R f = prefactor;
        for (int j = b; j < e; ++j) {
          int li_j = m.linear_index(th_, j);
          for (int i = b_in[r_]; i < e_in[r_]; ++i) {
            int li = li_j + i * m.stride(r_);
            output[li] += f * (E[comp][li] - E[comp][li - m.stride(th_)]);
          }
        }
      }

      if(hi and comp == phi_) {
        R f = prefactor / 3.0;
        for (int j = e; j < e+bdry; ++j) {
          int li_j = m.linear_index(th_, j-1);
          for (int i = b_in[r_]; i < e_in[r_]; ++i) {
            int li = li_j + i * m.stride(r_);
            auto x = E[comp][li] + E[comp][li];
            x += x;
            x += x;
            x += E[comp][li];
            output[li+m.stride(th_)] -= f * (x - E[comp][li - m.stride(th_)]);
          }
        }
      }

    }

  };

  template < typename R, int D >
  struct Diff_in_CurlB {
    const int mode; // 0 : d_r 1 : d_th
    const std::optional<int>& lo;
    const std::optional<int>& hi;

    const _helper<R,D>& h;

    const R dlnr;
    const R dth;

    using proxy_t =
      std::tuple<const int, const std::optional<int>&, const std::optional<int>&, const _helper<R,D>& , const Field<R, 3, D> &, const int, R, const int (&)[D], const int (&)[D]>;
    static constexpr int Prefactor = 6;

    proxy_t operator()(const Field<R, 3, D> &B, const int comp, R prefactor,
                       const int (&b)[D], const int (&e)[D]) {
      assert(comp != mode);
      prefactor /= (mode == r_ ? dlnr : dth);
      return std::forward_as_tuple(mode,lo,hi,h,B,comp, prefactor,b,e);
    }

    // \mathscr{h}_r^{-1} = \mathscr{h}_\theta^{-1} = 1 / (e^r \sin \theta),
    // \mathscr{h}_\phi^{-1} = \sin \theta / e^r
    // static constexpr bool o_B_off_diag = false;

    static void d_r_dB(Component<R,D,false>&& output, proxy_t&& diff) { // TODOL sematics
      const auto&&[mode,lo,hi,h,B,comp,prefactor,b_in,e_in] = std::move(diff);

      const auto& m = B.mesh();
      int b = b_in[r_];
      int e = e_in[r_];

      if(!hi) e--;

      R f = prefactor * h.estag() ;
      R g = h.estag() * h.estag();
      for (int j = b_in[th_]; j < e_in[th_]; ++j) {
        int li_j = m.linear_index(th_, j);
        for (int i = b; i < e; ++i) {
          int li = li_j + i * m.stride(r_);
          output[li] += f * h.fr(i) * ( B[comp][li + m.stride(r_)]*g - B[comp][li]);
        }
      }
      // NOTE turns out the forms of dr_dB_th and dr_dB_phi with scales are the same.

      // NOTE at outer boundary, the r derivative of B_\theta is half cell beyond axis, and r derivative of B_\phi is zero. Neither case needs explicit treatment
    }


    static void d_th_dB(Component<R,D,false>&& output, proxy_t&& diff) { // TODOL sematics
      const auto&&[mode, lo, hi,h,B,comp,prefactor,b_in,e_in] = std::move(diff);

      const auto& m = B.mesh();
      int b = b_in[th_];
      int e = e_in[th_] - 1;

      if ( comp == r_ ) {
        if (lo) {
          b = bdry;
          R f = prefactor / 3.0;
          f *= h.estag() * h.estag();
          for ( int j = 0; j < b; ++j ) {
            int li_j = m.linear_index(th_, j+1);
            R sr = h.sin(j + 2) / h.sin(j + 1);
            R fsr = f * h.sin(j, true) / h.sin(j + 2);
            for (int i = b_in[r_]; i < e_in[r_]; ++i) {
              int li = li_j + i * m.stride(r_);
              output[li-m.stride(th_)] += fsr * h.fr(i) * ( B[comp][li + m.stride(th_)] - B[comp][li] * sr );
            }
          }
        }
        {
          R f = prefactor;
          f *= h.estag() * h.estag();
          for (int j = b; j < e; ++j) {
            int li_j = m.linear_index(th_, j);
            R sr = h.sin(j + 1) / h.sin(j);
            R fsr = f * h.sin(j, true) / h.sin(j + 1);
            for (int i = b_in[r_]; i < e_in[r_]; ++i) {
              int li = li_j + i * m.stride(r_);
              output[li] += fsr * h.fr(i) * (B[r_][li + m.stride(th_)] - B[r_][li] * sr) ;
            }
          }
        }
        if (hi) {
          R f = prefactor / 3.0;
          f *= h.estag() * h.estag();
          for ( int j = e; j < e+bdry; ++j ) {
            int li_j = m.linear_index(th_, j);
            R sr = h.sin(j - 1) / h.sin(j);
            R fsr = f * h.sin(j, true) / h.sin(j - 1);
            for (int i = b_in[r_]; i < e_in[r_]; ++i) {
              int li = li_j + i * m.stride(r_);
              output[li] -= fsr * h.fr(i) * ( B[comp][li - m.stride(th_)] - B[comp][li]*sr );
            }
          }
        }
      }
      else { // comp = _phi
        // FIXME this part is not tested yet
        R f = prefactor;
        for (int j = b; j < e; ++j) {
          int li_j = m.linear_index(th_, j);
          R sinj0_in = h.sin(j);
          R sinj1_in = h.sin(j + 1);
          R fsin_inv = f / h.sin(j, true);
          for (int i = b_in[r_]; i < e_in[r_]; ++i) {
            int li = li_j + i * m.stride(r_);
            // FIXME which sin is nonzero?
            output[li] += fsin_inv * h.fr(i) *
              (B[phi_][li + m.stride(th_)] * sinj1_in -
               B[phi_][li] * sinj0_in);
          }
        }
      }
    }
  };

}

template <typename R, int D>
void operator+=( field::Component<R, D, false> &&output, typename field::Diff_in_CurlE<R,D>::proxy_t&& diff ) {
  using namespace field;
  const int mode = std::get<0>(diff);
  assert(mode < 2);
  if (mode == r_)
    Diff_in_CurlE<R,D>::d_r_dE(std::move(output), std::move(diff));
  else
    Diff_in_CurlE<R, D>::d_th_dE(std::move(output), std::move(diff));
}

template <typename R, int D>
void operator+=( field::Component<R, D, false> &&output, typename field::Diff_in_CurlB<R,D>::proxy_t&& diff ) {
  using namespace field;
  const int mode = std::get<0>(diff);
  assert(mode < 2);
  if (mode == r_)
    Diff_in_CurlB<R,D>::d_r_dB(std::move(output), std::move(diff));
  else
    Diff_in_CurlB<R, D>::d_th_dB(std::move(output), std::move(diff));
}

template <typename R, int D>
void operator-=(field::Component<R, D, false> &&output,
                typename field::Diff_in_CurlE<R, D>::proxy_t &&diff) {
  using namespace field;
  std::get<Diff_in_CurlE<R,D>::Prefactor>(diff) *= -1;
  std::move(output) += std::move(diff);
}

template <typename R, int D>
void operator-=(field::Component<R, D, false> &&output,
                typename field::Diff_in_CurlB<R, D>::proxy_t &&diff) {
  using namespace field;
  std::get<Diff_in_CurlB<R, D>::Prefactor>(diff) *= -1;
  std::move(output) += std::move(diff);
}

namespace field {
  template <typename R, int DGrid>
  class ortho_coord_transform_t {
  private:
    inline void _impl_no_sin( Component<R,DGrid,false> F, const int b_r, const int e_r, const int b_th, const int e_th, const _helper<R,DGrid>& h ) const { // TODOL semantics
      const auto& m = F.mesh();
      for (int j = b_th; j < e_th; ++j) {
        int li_j = m.linear_index(th_, j);
        for (int i = b_r; i < e_r; ++i) {
          int li = li_j + (i) * m.stride(r_);
          F[li] *= h.fr(i);
        }
      }
    }

    inline void _impl_sin( Component<R,DGrid,false> F, const int b_r, const int e_r, const int b_th, const int e_th, const _helper<R,DGrid>& h, bool ofs ) const { // TODOL semantics
      const auto& m = F.mesh();
      for (int j = b_th; j < e_th; ++j) {
        int li_j = m.linear_index(th_, j);
        R f = h.sin(j,ofs);
        for (int i = b_r; i < e_r; ++i) {
          int li = li_j + i * m.stride(r_);
          F[li] *= f * h.fr(i);
        }
      }
    }

    inline void _impl_sin_recip( Component<R,DGrid,false> F, const int b_r, const int e_r, const int b_th, const int e_th, const _helper<R,DGrid>& h) const { // TODOL semantics
      const auto& m = F. mesh();
      for (int j = b_th; j < e_th; ++j) {
        int li_j = m.linear_index(th_, j);
        R f = 1.0 / h.sin(j,true);
        for (int i = b_r; i < e_r; ++i) {
          int li = li_j + i * m.stride(r_);
          F[li] *= f * h.fr(i);
        }
      }
    }

    inline void _impl_c2o_B0( Component<R,DGrid,false> B0, const int b_r, const int e_r, const int b_th, const int e_th, const _helper<R,DGrid>& h) const { // TODOL semantics
      const auto& m = B0.mesh();
      R y = h.estag() * h.estag();
      for (int j = b_th; j < e_th; ++j) {
        int li_j = m.linear_index(th_, j);
        R f = y / h.sin(j);
        for (int i = b_r; i < e_r; ++i) {
          int li = li_j + i * m.stride(r_);
          B0[li] *= f * h.fr(i);
        }
      }
    }

  public:
    inline void ortho2coord( const std::string& name, Component<R,DGrid,false> F, const int b_r, const int e_r, const int b_th, const int e_th, const _helper<R,DGrid>& h ) const { // TODOL semantics
      if ( name == "E0" or name == "E1" or name == "B2")
        _impl_no_sin(F, b_r,e_r,b_th,e_th,h);
      else
        _impl_sin(F, b_r,e_r,b_th,e_th,h, (name != "B0"));
    }

    inline void coord2ortho( const std::string& name, Component<R,DGrid,false> F, const int b_r, const int e_r, const int b_th, const int e_th, const _helper<R,DGrid>& h ) const { // TODOL semantics
      if ( name == "E0" or name == "E1" or name == "B2")
        _impl_no_sin(F, b_r,e_r,b_th,e_th,h);
      else if ( name != "B0" )
        _impl_sin_recip(F, b_r,e_r,b_th,e_th,h);
      else
        _impl_c2o_B0(F, b_r,e_r,b_th,e_th,h );
    }
  };

  template <typename R, int D>
  inline auto check_boundaries(const apt::Grid<R,D>& grid, R lnr_surf, R lnr_out) {
    std::optional<int> surf, outer, lo_axis, hi_axis;
    if (grid[0].lower() < lnr_surf + grid[0].delta() and grid[0].upper() > lnr_surf)
      surf.emplace(grid[0].csba(lnr_surf));
    if (grid[0].lower() < lnr_out and grid[0].upper() > lnr_out - grid[0].delta())
      outer.emplace(grid[0].csba(lnr_out));
    if (grid[1].lower() < grid[1].delta())
      lo_axis.emplace(0);
    if (grid[1].upper() > std::acos(-1.0) - grid[1].delta())
      hi_axis.emplace(grid[1].dim());
    return std::make_tuple(surf,outer,lo_axis,hi_axis);
  }

  template <typename R, int DGrid, typename RJ>
  void LogSphericalSolver<R, DGrid, RJ>
  ::operator()(Field<R, 3, DGrid> &E,
               Field<R, 3, DGrid> &B,
               Field<RJ, 3, DGrid> &J,
               const apt::Grid<R, DGrid> &grid,
               const mpi::CartComm &comm,
               int timestep, R dt) const {
    // FIXME
    // auto check_invalid = [](const auto &f, auto name) {
    //   for (int c = 0; c < 3; ++c) {
    //     int res = 0;
    //     for (const auto &x : f[c].data())
    //       res += (std::isnan(x) or std::abs(x) > 1e10);
    //     if (res > 0) {
    //       std::cout << "FOUND INVALID! " << name << c << std::endl;
    //     }
    //   }
    // };

    static_assert(DGrid == 2);
    // NOTE E,B,J are assumed to have been synced, and have same mesh

    const R alpha = std::min<R>(_alpha, 1.0);
    const R beta = 1.0 - alpha;

    const auto[surf, outer, lo_axis, hi_axis] = check_boundaries(grid,_surface,_outer);

    const auto& mesh = E.mesh();
    int b[DGrid] = { (surf ? *surf : mesh.range(r_).far_begin()),
                     (lo_axis ? *lo_axis : mesh.range(th_).far_begin() )   };
    int e[DGrid] = { (outer ? *outer : mesh.range(r_).far_end()),
                     (hi_axis ? *hi_axis : mesh.range(th_).far_end())   };

    if (!hopt<R,DGrid>)
      hopt<R,DGrid>.emplace(grid, b[r_] - bool(surf), e[r_], b[th_], e[th_]);

    const auto &h = *(hopt<R,DGrid>);
    const ortho_coord_transform_t<R,DGrid> tr;
    {
      for ( int i = b[r_]-bool(surf); i < e[r_]; ++i )
        h.fr(i) = std::exp(grid[r_].absc(i));
      // FIXME turn back tr on
      // tr.ortho2coord("E0", E[0], b[r_], e[r_], b[th_], e[th_], h);
      for ( auto& x : h.fr() ) x /= h.estag();
      // tr.ortho2coord("E1", E[1], b[r_], e[r_], b[th_], e[th_], h);
      // tr.ortho2coord("E2", E[2], b[r_], e[r_], b[th_], e[th_], h);
      for (auto &x : h.fr() ) x *= x;
      // tr.ortho2coord("B0", B[0], b[r_], e[r_], b[th_], e[th_], h); // FIXME check boundary
      R y = h.estag() * h.estag();
      for (auto &x : h.fr()) x *= y;
      // tr.ortho2coord("B1", B[1], b[r_], e[r_], b[th_], e[th_], h);
      // tr.ortho2coord("B2", B[2], b[r_], e[r_], b[th_], e[th_], h);
      // NOTE assert: now h.fr stores (r_i)^2 at integer cells
    }

    // ortho2coord(E, B, grid, b[r_]-bool(surf), e[r_], b[th_], e[th_], h);
    for (auto &x : h.fr()) x = 1.0 / x;

    const auto dlnr = grid[r_].delta();
    const auto dth = grid[th_].delta();

    Diff_in_CurlE<R, DGrid> diff_r{r_, surf, outer, dlnr, dth};
    Diff_in_CurlE<R, DGrid> diff_th{th_, lo_axis, hi_axis, dlnr, dth};
    Diff_in_CurlB<R, DGrid> diffh_r{r_, surf, outer, h, dlnr, dth};
    Diff_in_CurlB<R, DGrid> diffh_th{th_, lo_axis, hi_axis, h, dlnr, dth};

    auto after_curlE =
      [&surf,&lo_axis]( auto& b, auto& e ) {
        if (!surf) b[r_]++;
        if (!lo_axis) b[th_]++;
      };
    auto after_curlB =
      [&outer, &hi_axis](auto &b, auto &e) {
        if (!outer) e[r_]--;
        if (!hi_axis) e[th_]--;
      };

    auto mul_absorb =
      []( auto& output, const auto& input, auto factor, const auto&b, const auto& e ){
        const auto& m = input.mesh();
        for ( int c = 0; c < 3; ++c ) {
          for (int j = b[th_]; j < e[th_]; ++j) {
            int li_j = m.linear_index(th_,j);
            for (int i = b[r_]; i < e[r_]; ++i) {
              int li = li_j + i * m.stride(r_);
              output[c][li] += factor * input[c][li];
            }
          }
        }
      };

    auto absorb_curl
      = [&diff_r, &diff_th,&after_curlE](auto &F, const auto &E, auto factor, auto &b, auto &e) {
      F[r_] += diff_th(E, phi_, factor, b, e);
      F[th_] -= diff_r(E, phi_, factor, b, e);
      F[phi_] += diff_r(E, th_, factor, b, e);
      F[phi_] -= diff_th(E, r_, factor, b, e);
      after_curlE(b, e);
    };

    auto absorb_curlh
      = [&diffh_r, &diffh_th,&after_curlB](auto &F, const auto &B, auto factor, auto &b, auto &e) {
      F[r_] += diffh_th(B, phi_, factor, b, e);
      F[th_] -= diffh_r(B, phi_, factor, b, e);
      F[phi_] += diffh_r(B, th_, factor, b, e);
      F[phi_] -= diffh_th(B, r_, factor, b, e);
      after_curlB(b, e);
    };

    const bool beta_neq_0 = (beta > 1e-8);

    Field<R, 3, DGrid> C (mesh);
    {
      C[r_] += diff_th(E,phi_,1.0,b,e);
      C[th_] -= diff_r(E,phi_,1.0,b,e,true);
      C[phi_] += diff_r(E,th_,1.0,b,e,true);
      C[phi_] -= diff_th(E,r_,1.0,b,e);
    }

    if (beta_neq_0) {
      after_curlE(b, e);
      mul_absorb(B, C, -alpha * beta * dt, b, e);
    }
    for ( int c = 0; c < 3; ++c ) {
      auto f = -4.0*std::acos(-1.0)*dt;
      for ( auto& x : J[c].data() )
        x *= f;
    }
    absorb_curlh ( J, B, dt, b, e );
    mul_absorb(E, J, 1.0, b, e );
    if (beta_neq_0)
      mul_absorb(B, C, -beta*beta*dt, b,e);

    for (int i = std::max(_op_inv_precision, 1); i != 0; --i) {
      absorb_curl(C,J,1.0,b,e);
      if ( i != 1 ) {
        J.reset();
        absorb_curlh(J, C, -(alpha * alpha * dt * dt), b, e);
        mul_absorb(E, J, 1.0, b, e);
        C.reset(); // crucial to place it here
      } else {
        // last iteration
        absorb_curlh(E, C, -(alpha * alpha * dt * dt), b, e);
      }
    }

    absorb_curl(B,E,-(alpha * dt), b,e);
    // NOTE sync E B is done outside

    // {
    //   //NOTE: assert: h.fr contains 1 / (r_i.0)^2
    //   tr.coord2ortho("B1", B[1], b[r_], e[r_], b[th_], e[th_], h);
    //   tr.coord2ortho("B2", B[2], b[r_], e[r_], b[th_], e[th_], h);
    //   for ( auto& x : h.fr() ) x = std::sqrt(x);
    //   tr.coord2ortho("E0", E[0], b[r_], e[r_], b[th_], e[th_], h);
    //   for (auto &x : h.fr() ) x *= h.estag();
    //   tr.coord2ortho("E1", E[1], b[r_], e[r_], b[th_], e[th_], h);
    //   tr.coord2ortho("E2", E[2], b[r_], e[r_], b[th_], e[th_], h);
    //   tr.coord2ortho("B0", B[0], b[r_], e[r_], b[th_], e[th_], h); // FIXME check boundary
    //   // NOTE assert: now h.fr stores (r_i)^2 at integer cells
    // }

    // coord2ortho(B, E, grid, b[r_] - bool(surf), e[r_], b[th_], e[th_], h); // FIXME double using the shrunken b,e is OK. OK for bulk // patches, how about for surface?
  }
}
