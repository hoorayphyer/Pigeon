#include "field/log_spherical_solver/updater.hpp"
// #include "field/sync.hpp" // FIXME not in use right now
#include "field/yee.hpp"
#include <optional>
#include <cassert>
#include <tuple>

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
  private:
    static void _diff_bulk( Component<R,D,false>& output,
                            const Component<R,D>& input,
                            const R f, const int dir,
                            const int b_r, const int e_r,
                            const int b_th, const int e_th
                            ) {
      const auto& m = input.mesh();
      const int s = m.stride(dir);
      for (int j = b_th; j < e_th; ++j) {
        int li_j = m.linear_index(th_, j);
        for (int i = b_r; i < e_r; ++i) {
          int li = li_j + i * m.stride(r_);
          output[li] += f * (input[li] - input[li - s]);
        }
      }
    }

  public:
    const int mode; // 0 : d_r 1 : d_th
    const std::optional<int>& lo;
    const std::optional<int>& hi;
    const R delta;

    using proxy_t =
        std::tuple<const int, const std::optional<int>&, const std::optional<int>&,
                   const Field<R,3,D> &, const int, R,
                   const int (&)[D], const int (&)[D], const bool >;
    static constexpr int Prefactor = 5;

    proxy_t operator()(const Field<R, 3, D> &E, const int comp, R prefactor,
                       const int (&b)[D], const int (&e)[D],
                       const bool continuous = false) {
      assert(comp != mode);
      prefactor /= delta;
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

        const auto& s = m.stride(r_);
        for ( int j = b_in[th_]; j < e_in[th_]; ++j ) {
          int li_j = m.linear_index(th_, j);
          for ( int i = *lo; i < b; ++i ) {
            int li = li_j + (i+1) * s;
            // The following is O(dr^2) precise.
            auto x = E[comp][li] - E[comp][li-s];
            output[li-s] += f * ( x + x + E[comp][li] - E[comp][li+s] );
          }
        }
      }

      _diff_bulk(output, E[comp], prefactor, r_, b, e, b_in[th_], e_in[th_] );
      // NOTE at outer boundary, the r derivative on E_\theta and E_\phi is taken to vanish.
    }

    static void d_th_dE_r(Component<R,D,false>&& output, proxy_t&& diff) { // TODOL sematics
      const auto &&[mode, lo, hi, E, comp, prefactor, b_in, e_in, continuous] = std::move(diff);
      int b = b_in[th_];
      int e = e_in[th_] - bdry*bool(hi);

      if(!lo) b++;
      else {
        b = bdry; // NOTE true for both comp == phi_ and comp == r_
        // NOTE if comp == r_, by symmetry output += 0, which doesn't need any action.
      }

      _diff_bulk(output, E[comp], prefactor, th_, b_in[r_], e_in[r_],b, e );
    }

    static void d_th_dE_phi(Component<R,D,false>&& output, proxy_t&& diff) { // TODOL sematics
      const auto &&[mode, lo, hi, E, comp, prefactor, b_in, e_in, continuous] = std::move(diff);

      int b = b_in[th_];
      int e = e_in[th_] - bdry*bool(hi);

      auto diff_bdry
        = [&, &m=E.mesh()]( const R f, const int j, const bool hi ) {
          static_assert(bdry == 1);
          const int li_j = m.linear_index(th_, j);
          const int y = hi ? -1 : 1;
          for (int i = b_in[r_]; i < e_in[r_]; ++i) {
            int li = li_j + i * m.stride(r_);
            auto x = E[phi_][li] + E[phi_][li];
            x += x;
            x += x;
            x += E[phi_][li];
            output[li + (hi)*m.stride(th_)] += f * ( x - E[phi_][li + y*m.stride(th_)] );
          }
        };

      if(!lo) b++;
      else {
        b = bdry; // NOTE true for both comp == phi_ and comp == r_
        diff_bdry(prefactor/3.0, 0, false);
      }

      _diff_bulk(output, E[comp], prefactor, th_, b_in[r_], e_in[r_],b, e );

      if(hi)
        diff_bdry(-prefactor/3.0, e-1, true);
    }

  };

  template < typename R, int D >
  struct Diff_in_CurlB {
    const int mode; // 0 : d_r 1 : d_th
    const std::optional<int>& lo;
    const std::optional<int>& hi;
    const _helper<R,D>& h;
    const R delta;

    using proxy_t =
      std::tuple<const int, const std::optional<int>&, const std::optional<int>&, const _helper<R,D>& , const Field<R, 3, D> &, const int, R, const int (&)[D], const int (&)[D]>;
    static constexpr int Prefactor = 6;

    proxy_t operator()(const Field<R, 3, D> &B, const int comp, R prefactor,
                       const int (&b)[D], const int (&e)[D]) {
      assert(comp != mode);
      prefactor /= delta;
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


    static void d_th_dB_r(Component<R,D,false>&& output, proxy_t&& diff) { // TODOL sematics
      const auto&&[mode, lo, hi,h,B,comp,prefactor,b_in,e_in] = std::move(diff);

      const auto& m = B.mesh();
      int b = b_in[th_];
      int e = e_in[th_] - 1 - bdry*bool(hi);

      auto diff_bdry
        = [&, &m=B.mesh()]( const R f, const int j, const bool lo ) {
          static_assert(bdry == 1);
          const int li_j = m.linear_index(th_, j);
          const int y = lo ? 1 : -1;
          const R sr = h.sin(j + y) / h.sin(j);
          const R fsr = f * h.sin(j, true) / h.sin(j + y);

          for (int i = b_in[r_]; i < e_in[r_]; ++i) {
            int li = li_j + i * m.stride(r_);
            output[li-(lo)*m.stride(th_)] += fsr * h.fr(i) * ( B[comp][li + y * m.stride(th_)] - B[comp][li] * sr );
          }
        };

      if (lo) {
        b = bdry;
        diff_bdry( prefactor*h.estag()*h.estag()/3.0, 1, true );
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
      if (hi)
        diff_bdry(-prefactor * h.estag() * h.estag() / 3.0, e, false);
    }

    static void d_th_dB_phi(Component<R,D,false>&& output, proxy_t&& diff) { // TODOL sematics
      const auto&&[mode, lo, hi,h,B,comp,prefactor,b_in,e_in] = std::move(diff);

      const auto& m = B.mesh();
      const int b = b_in[th_];
      const int e = e_in[th_] - 1;

      const R f = prefactor;
      for (int j = b; j < e; ++j) {
        const int li_j = m.linear_index(th_, j);
        const R sinj0_in = h.sin(j);
        const R sinj1_in = h.sin(j + 1);
        const R fsin_inv = f / h.sin(j, true);
        for (int i = b_in[r_]; i < e_in[r_]; ++i) {
          int li = li_j + i * m.stride(r_);
          output[li] += fsin_inv * h.fr(i) *
            (B[phi_][li + m.stride(th_)] * sinj1_in -
             B[phi_][li] * sinj0_in);
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
  const int comp = std::get<5>(diff);
  if (mode == r_)
    Diff_in_CurlE<R,D>::d_r_dE(std::move(output), std::move(diff));
  else if (comp == r_)
    Diff_in_CurlE<R, D>::d_th_dE_r(std::move(output), std::move(diff));
  else
    Diff_in_CurlE<R, D>::d_th_dE_phi(std::move(output), std::move(diff));
}

template <typename R, int D>
void operator+=( field::Component<R, D, false> &&output, typename field::Diff_in_CurlB<R,D>::proxy_t&& diff ) {
  using namespace field;
  const int mode = std::get<0>(diff);
  assert(mode < 2);
  const int comp = std::get<5>(diff);
  assert(comp == r_ or comp == phi_);
  if (mode == r_)
    Diff_in_CurlB<R,D>::d_r_dB(std::move(output), std::move(diff));
  else if ( comp == r_ )
    Diff_in_CurlB<R, D>::d_th_dB_r(std::move(output), std::move(diff));
  else
    Diff_in_CurlB<R, D>::d_th_dB_phi(std::move(output), std::move(diff));
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
          F[li_j + i*m.stride(r_)] *= h.fr(i);
        }
      }
    }

    inline void _impl_sin( Component<R,DGrid,false> F, const int b_r, const int e_r, const int b_th, const int e_th, const _helper<R,DGrid>& h, const bool ofs_th ) const { // TODOL semantics
      const auto& m = F.mesh();
      for (int j = b_th; j < e_th; ++j) {
        int li_j = m.linear_index(th_, j);
        R f = h.sin(j,ofs_th);
        for (int i = b_r; i < e_r; ++i) {
          F[li_j + i * m.stride(r_)] *= f * h.fr(i);
        }
      }
    }

    inline void _impl_sin_recip( Component<R,DGrid,false> F, const int b_r, const int e_r, const int b_th, const int e_th, const _helper<R,DGrid>& h) const { // TODOL semantics
      const auto& m = F. mesh();
      for (int j = b_th; j < e_th; ++j) {
        int li_j = m.linear_index(th_, j);
        R f = 1.0 / h.sin(j,true);
        for (int i = b_r; i < e_r; ++i) {
          F[li_j + i * m.stride(r_)] *= f * h.fr(i);
        }
      }
    }

    inline void _impl_c2o_B0( Component<R,DGrid,false> B0, const int b_r, const int e_r, const int b_th, const int e_th, const _helper<R,DGrid>& h, const bool lo_axis, const bool hi_axis) const { // TODOL semantics
      const auto& m = B0.mesh();
      R y = h.estag() * h.estag();
      auto f_axis = [&]( int j, bool lo ) {
        int sign = lo ? 1 : -1;
        j += sign;
        int li_j = m.linear_index(th_, j);
        int s = sign * m.stride(th_);
        const R f = y / (6.0 * std::asin(h.sin(j)) ); // asin(sin(j)) is a trick to get delta theta
        for (int i = b_r; i < e_r; ++i) {
          int li = li_j + i * m.stride(r_);
          auto x = B0[li] + B0[li];
          x += x;
          x += x;
          B0[li - s] = f * h.fr(i) * (x - B0[li + s]);
        }
      };

      if ( lo_axis ) f_axis(b_th, true);
      if (hi_axis) f_axis(e_th-1, false); // NOTE it uses coordinate B0, so must go before the bulk

      for (int j = b_th+lo_axis; j < e_th-hi_axis; ++j) {
        int li_j = m.linear_index(th_, j);
        R f = y / h.sin(j);
        for (int i = b_r; i < e_r; ++i) {
          B0[li_j + i * m.stride(r_)] *= f * h.fr(i);
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

    inline void coord2ortho( const std::string& name, Component<R,DGrid,false> F, const int b_r, const int e_r, const int b_th, const int e_th, const _helper<R,DGrid>& h, const bool lo_axis=false, const bool hi_axis=false ) const { // TODOL semantics
      if ( name == "E0" or name == "E1" or name == "B2")
        _impl_no_sin(F, b_r,e_r,b_th,e_th,h);
      else if ( name != "B0" )
        _impl_sin_recip(F, b_r,e_r,b_th,e_th,h);
      else
        _impl_c2o_B0(F, b_r,e_r,b_th,e_th,h,lo_axis,hi_axis );
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
      hi_axis.emplace(grid[1].dim() + 1);
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
    // FIXME not using range

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
    // NOTE the hi_axis is included

    if (!hopt<R,DGrid>)
      hopt<R,DGrid>.emplace(grid, b[r_] - bool(surf), e[r_], b[th_], e[th_]);

    const auto &h = *(hopt<R,DGrid>);
    const ortho_coord_transform_t<R,DGrid> tr;
    {
      for ( int i = b[r_]-bool(surf); i < e[r_]; ++i )
        h.fr(i) = std::exp(grid[r_].absc(i));
      tr.ortho2coord("E0", E[0], b[r_], e[r_], b[th_], e[th_]-bool(hi_axis), h);
      for ( auto& x : h.fr() ) x /= h.estag();
      tr.ortho2coord("E1", E[1], b[r_]-bool(surf), e[r_], b[th_], e[th_], h);
      tr.ortho2coord("E2", E[2], b[r_]-bool(surf), e[r_], b[th_], e[th_]-bool(hi_axis), h);
      for (auto &x : h.fr() ) x *= x;
      tr.ortho2coord("B0", B[0], b[r_], e[r_], b[th_], e[th_], h);
      R y = h.estag() * h.estag();
      for (auto &x : h.fr()) x *= y;
      tr.ortho2coord("B1", B[1], b[r_], e[r_], b[th_], e[th_]-bool(hi_axis), h);
      tr.ortho2coord("B2", B[2], b[r_], e[r_], b[th_], e[th_], h);
      // NOTE assert: now h.fr stores (r_i.0)^2
    }

    for (auto &x : h.fr()) x = 1.0 / x;

    const auto dlnr = grid[r_].delta();
    const auto dth = grid[th_].delta();

    Diff_in_CurlE<R, DGrid> diff_r{r_, surf, outer, dlnr};
    Diff_in_CurlE<R, DGrid> diff_th{th_, lo_axis, hi_axis, dth};
    Diff_in_CurlB<R, DGrid> diffh_r{r_, surf, outer, h, dlnr};
    Diff_in_CurlB<R, DGrid> diffh_th{th_, lo_axis, hi_axis, h, dth};

    auto after_curlE =
      [&surf,&lo_axis]( auto& b) {
        if (!surf) b[r_]++;
        if (!lo_axis) b[th_]++;
      };
    auto after_curlB =
      [&outer, &hi_axis](auto &e) {
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
      = [&diff_r, &diff_th,&after_curlE,hi=bool(hi_axis)](auto &F, const auto &E, auto factor, auto &b, auto &e, bool continuous = false) {
      F[r_] += diff_th(E, phi_, factor, b, e);
      e[th_] -= hi;
      F[th_] -= diff_r(E, phi_, factor, b, e, continuous);
      e[th_] += hi;
      F[phi_] += diff_r(E, th_, factor, b, e, continuous);
      F[phi_] -= diff_th(E, r_, factor, b, e);
      after_curlE(b);
    };

    auto absorb_curlh
      = [&diffh_r, &diffh_th,&after_curlB,hi=bool(hi_axis)](auto &F, const auto &B, auto factor, auto &b, auto &e) {
      F[r_] += diffh_th(B, phi_, factor, b, e);
      F[th_] -= diffh_r(B, phi_, factor, b, e);
      e[th_] -= hi;
      F[phi_] += diffh_r(B, th_, factor, b, e);
      e[th_] += hi;
      F[phi_] -= diffh_th(B, r_, factor, b, e);
      after_curlB(e);
    };

    const bool beta_neq_0 = (beta > 1e-8);

    for (int c = 0; c < 3; ++c) {
      auto f = -_fourpi * dt;
      for (auto &x : J[c].data())
        x *= f;
    }

    Field<R, 3, DGrid> C (mesh);
    absorb_curl(C,E,1.0, b,e, true); // Eq(*)

    if (beta_neq_0) {
      mul_absorb(B, C, -alpha * beta * dt, b, e);
    } else { // undo after_curlE in Eq(*) for the sake of Eq(**)
      if (!surf) b[r_]--;
      if (!lo_axis) b[th_]--;
    }
    absorb_curlh ( J, B, dt, b, e ); // Eq(**)
    mul_absorb(E, J, 1.0, b, e );
    if (beta_neq_0) mul_absorb(B, C, -beta*beta*dt, b,e);

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

    {
      // NOTE: assert: h.fr contains 1 / (r_i.0)^2
      tr.coord2ortho("B0", B[0], b[r_], e[r_], b[th_], e[th_], h, bool(lo_axis), bool(hi_axis));
      tr.coord2ortho("B1", B[1], b[r_], e[r_], b[th_], e[th_], h);
      tr.coord2ortho("B2", B[2], b[r_], e[r_], b[th_], e[th_]-bool(hi_axis), h);
      for ( auto& x : h.fr() ) x = std::sqrt(x);
      tr.coord2ortho("E0", E[0], b[r_], e[r_], b[th_], e[th_]-bool(hi_axis), h);
      for (auto &x : h.fr() ) x *= h.estag();
      tr.coord2ortho("E1", E[1], b[r_]-bool(surf), e[r_], b[th_], e[th_], h);
      tr.coord2ortho("E2", E[2], b[r_]-bool(surf), e[r_], b[th_], e[th_]-bool(hi_axis), h);
      // NOTE assert: now h.fr stores (r_i)^2 at integer cells
    }

  }
}
