#include "field/log_spherical_solver/updater.hpp"
#include <optional>
#include <cassert>
#include <tuple>

namespace field {
  constexpr int r_ = 0, th_ = 1, phi_ = 2;
  constexpr int bdry = 1; // set by 2nd order scheme

  template <typename R, int D>
  class _helper {
  private:
    std::vector<R> _fr{};
    std::vector<R> _sin_th{};
    int _b_r;
    int _b_th;
    R _estag; // exp( 0.5 * dlnr )

  public:
    _helper() = default;
    _helper(const apt::Grid<R, D> &grid, int b_r, int e_r, int b_th, int e_th)
        : _b_r(b_r), _b_th(b_th), _estag(std::exp(grid[r_].delta()*0.5)) {
      _fr.resize(e_r - b_r);
      _fr.shrink_to_fit();
      _sin_th.resize(2 * (e_th - b_th));
      _sin_th.shrink_to_fit();

      for (int i = 0; i < e_r - b_r; ++i)
        _fr[i] = std::exp(grid[r_].absc(b_r+i));

      for (int j = 0; j < e_th - b_th; ++j) {
        _sin_th[2 * j] = std::sin(grid[th_].absc(b_th + j));
        _sin_th[2 * j + 1] = std::sin(grid[th_].absc(b_th + j, 0.5));
      }
    }

    const std::vector<R> &fr() const { return _fr; }
    std::vector<R> &fr() { return _fr; }
    const auto &fr(int i) const { return _fr[i - _b_r]; }
    auto &fr(int i) { return _fr[i - _b_r]; }

    std::vector<R> &fsin() { return _sin_th; }

    const auto &fsin(int j, bool midway = false) const {
      return _sin_th[((j - _b_th) << 1) + midway];
    }
    const auto& estag() const {return _estag;}
  };

  template < typename R, int D >
  struct Diff_in_CurlB {
    const int mode; // 0 : d_r 1 : d_th
    const std::optional<int>& lo;
    const std::optional<int>& hi;
    const R delta;

    void d(Component<R,D,false>&& output, const Field<R, 3, D> &B, const int comp, R prefactor, const int (&b_in)[D], const int (&e_in)[D]) {// TODOL sematics
      prefactor /= delta;
      const auto& m = B.mesh();

      if (r_ == mode and lo) { // conductor surface
        const auto &s = m.stride(r_);
        for (int j = b_in[th_]; j < e_in[th_]; ++j) {
          int li = m.linear_index(th_, j) + b_in[r_] * s;
          auto x = B[comp][li + s] - B[comp][li];
          output[li - s] +=
            prefactor * (x + x + B[comp][li + s] - B[comp][li + (s << 1)]);
        }
      }

      { // bulk
        const int s = m.stride(mode);
        for (int j = b_in[th_]; j < e_in[th_] - (th_==mode); ++j) {
          int li_j = m.linear_index(th_, j);
          for (int i = b_in[r_]; i < e_in[r_] - (r_==mode); ++i) {
            int li = li_j + i * m.stride(r_);
            output[li] += prefactor * (B[comp][li + s] - B[comp][li]);
          }
        }
      }

      // at r outer boundary, the r derivative is taken to vanish.
    }

  };

  template < typename R, int D >
  struct Diff_in_CurlE {
    const int mode; // 0 : d_r 1 : d_th
    const std::optional<int>& lo;
    const std::optional<int>& hi;
    const _helper<R,D>& h;
    const R delta;

    void d_r(Component<R,D,false>&& output, const Field<R, 3, D> &E, const int comp, R prefactor, const int (&b_in)[D], const int (&e_in)[D]) { // TODOL sematics
      prefactor /= delta;

      const auto& m = E.mesh();

      R f = prefactor / h.estag();
      R g = h.estag() * h.estag();
      for (int j = b_in[th_]; j < e_in[th_]; ++j) {
        int li_j = m.linear_index(th_, j);
        for (int i = b_in[r_] + 1; i < e_in[r_]; ++i) {
          int li = li_j + i * m.stride(r_);
          output[li] += f * h.fr(i) * ( E[comp][li] - E[comp][li - m.stride(r_)]*g );
        }
      }
    }

    void d_th(Component<R,D,false>&& output, const Field<R, 3, D> &E, const int comp, R prefactor, const int (&b_in)[D], const int (&e_in)[D]) { // TODOL sematics
      prefactor /= delta;

      const int b = b_in[th_];
      const int e = e_in[th_];

      const R f = (phi_ == comp) ? prefactor /( h.estag()* h.estag() ) : prefactor;

      auto diff_bdry
        = [&output, &e, &h=this->h, &E, &comp, b_r=b_in[r_], e_r=e_in[r_], &f]( const bool lo ) {
          const auto& m = E.mesh();

          const R ff = lo ? 4.0 * f : -4.0 * f;
          const int y = lo ? 0 : -1;

          const int li_j = m.linear_index(th_, lo ? 0 : e);
          for (int i = b_r; i < e_r; ++i) {
            int li = li_j + i * m.stride(r_);
            // NOTE order O(dth^2)
            output[li] += ff * h.fr(i) * E[comp][li + y * m.stride(th_)];
          }
        };

      if (phi_ == comp and lo) diff_bdry(true);

      const auto& m = E.mesh();
      const bool inv = (r_ == comp);
      for (int j = b+1; j < e; ++j) {
        int li_j = m.linear_index(th_, j);
        R fsr = f * h.fsin(j, !inv) / h.fsin(j, inv);
        R sr = h.fsin(j - (!inv), true) / h.fsin(j - inv, true);
        for (int i = b_in[r_]; i < e_in[r_]; ++i) {
          int li = li_j + i * m.stride(r_);
          output[li] += fsr * h.fr(i) * (E[comp][li] - E[comp][li - m.stride(th_)] * sr);
        }
      }

      if (phi_ == comp and hi) diff_bdry(false);
    }
  };

}

namespace field {
  template <typename R, int DGrid>
  void transform( const std::string& name, Component<R,DGrid,false> F, const int b_r, const int e_r, const int b_th, const int e_th, const _helper<R,DGrid>& h ) { // TODOL semantics
    auto impl = [&](const bool has_sin, const bool ofs_th) {
      const auto &m = F.mesh();
      for (int j = b_th; j < e_th; ++j) {
        int li_j = m.linear_index(th_, j);
        R f = has_sin ? h.fsin(j, ofs_th) : 1.0;
        for (int i = b_r; i < e_r; ++i) {
          F[li_j + i * m.stride(r_)] *= f * h.fr(i);
        }
      }
    };

    if ( name == "B0" or name == "B1" or name == "E2")
      impl(false, false);
    else
      impl(true, (name == "E0"));
  }

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
    static_assert(DGrid == 2);
    // NOTE E,B,J are assumed to have been synced, and have same mesh

    const R alpha = std::min<R>(_alpha, 1.0);

    const auto[surf, outer, lo_axis, hi_axis] = check_boundaries(grid,_surface,_outer);

    const auto& mesh = E.mesh(); // FIXME this is not using range
    int b[DGrid] = { (surf ? *surf : mesh.range(r_).far_begin()),
                     (lo_axis ? *lo_axis : mesh.range(th_).far_begin() )   };
    int e[DGrid] = { (outer ? *outer : mesh.range(r_).far_end()),
                     (hi_axis ? *hi_axis : mesh.range(th_).far_end())   };
    // NOTE b,e are defined as that of B_\phi

    const int e_r_max = outer ? mesh.range(r_).far_end() : e[r_];
    _helper<R, DGrid> h{grid, b[r_] - bool(surf), e_r_max, b[th_], e[th_]};

    auto tr = [&](R estag, const int e_r_max) {
      transform("B1", B[1], b[r_], e_r_max, b[th_], e[th_] - bool(hi_axis), h);
      transform("B2", B[2], b[r_], e_r_max, b[th_], e[th_], h);
      for (auto &x : h.fr()) x *= estag;
      transform("B0", B[0], b[r_] - bool(surf), e_r_max, b[th_], e[th_], h);

      for (auto &x : h.fr()) x *= x;
      transform("E1", E[1], b[r_] - bool(surf), e_r_max, b[th_], e[th_], h);
      transform("E2", E[2], b[r_] - bool(surf), e_r_max, b[th_], e[th_] - bool(hi_axis), h);
      estag *= estag;
      for (auto &x : h.fr()) x /= estag;
      transform("E0", E[0], b[r_], e_r_max, b[th_], e[th_] - bool(hi_axis), h);
    };
    tr(h.estag(), e_r_max);
    for (auto &x : h.fr()) x = 1.0 / x;
    // NOTE assert: now h.fr stores (r_i.0)^-2

    Diff_in_CurlB<R, DGrid> diff_r{r_, surf, outer, grid[r_].delta()};
    Diff_in_CurlB<R, DGrid> diff_th{th_, lo_axis, hi_axis, grid[th_].delta()};
    Diff_in_CurlE<R, DGrid> diffk_r{r_, surf, outer, h, grid[r_].delta()};
    Diff_in_CurlE<R, DGrid> diffk_th{th_, lo_axis, hi_axis, h, grid[th_].delta()};

    auto mul_absorb =
      []( auto& output, const auto& input, auto factor ){
        // assume input and output have same mesh
        const auto& m = input.mesh();
        for ( int c = 0; c < 3; ++c ) {
          for ( std::size_t i = 0; i < input[c].data().size(); ++i )
            output[c].data()[i] += factor * input[c].data()[i];
        }
      };

    auto absorb_curlk
      = [&diffk_r, &diffk_th,lo=bool(lo_axis),hi=bool(hi_axis),surf=bool(surf)](auto &F, const auto &E, const auto factor, auto &b, auto &e) {
      b[r_] -= surf;
      e[th_] -= hi;
      diffk_th.d_th(F[r_], E, phi_, factor, b, e);
      diffk_r.d_r(F[th_], E, phi_, -factor, b, e);
      e[th_] += hi;
      diffk_r.d_r(F[phi_], E, th_, factor, b, e); // NOTE d_r E_th is NOT enforced to be zero on higher axis
      b[r_] += surf;
      e[th_] -= hi;
      diffk_th.d_th(F[phi_], E, r_, -factor, b, e);
      e[th_] += hi;

      if (!surf) b[r_]++;
      if (!lo) b[th_]++;
    };

    auto absorb_curl
      = [&diff_r, &diff_th,hi=bool(hi_axis),surf=bool(surf),outer=bool(outer)](auto &F, const auto &B, const auto factor, auto &b, auto &e) {
      diff_th.d(F[r_], B, phi_, factor, b, e);
      diff_r.d(F[th_], B, phi_, -factor, b, e); // NOTE d_r B_phi is NOT enforced to be zero on higher axis
      e[th_] -= hi;
      diff_r.d(F[phi_], B, th_, factor, b, e);
      e[th_] += hi;
      b[r_] -=surf;
      diff_th.d(F[phi_], B, r_, -factor, b, e);
      b[r_] +=surf;

      if (!outer) e[r_]--;
      if (!hi) e[th_]--;
    };

    for (int c = 0; c < 3; ++c) {
      auto f = -_fourpi * dt;
      for (auto &x : J[c].data()) x *= f;
    }

    Field<R, 3, DGrid> C (mesh);
    absorb_curlk(C,E,1.0, b,e);

    mul_absorb(B, C, -alpha * dt);
    absorb_curl( J, B, dt, b, e );
    mul_absorb(B, C, (alpha-(1.0 - alpha)) * dt);

    mul_absorb(E, J, 1.0);
    for (int i = std::max(_op_inv_precision, 1); i != 1; --i) {
      C.reset();
      absorb_curlk(C,J,1.0,b,e);
      J.reset();
      absorb_curl(J, C, -(alpha * alpha * dt * dt), b, e);
      mul_absorb(E, J, 1.0);
    }
    C.reset();
    absorb_curlk(C,J,1.0,b,e);
    absorb_curl(E, C, -(alpha * alpha * dt * dt), b, e);

    absorb_curlk(B,E,-(alpha * dt), b,e);
    // NOTE sync E B is done outside

    {
      // NOTE: assert: h.fr contains 1 / (r_i.0)^2
      for (auto &x : h.fr()) x = std::sqrt(x);
      for (auto &x : h.fsin()) {
        x = (std::abs(x) < 0.5 * 0.5 * grid[th_].delta()) ? 0.0 : 1.0 / x;
      }
      tr(1/h.estag(), outer ? mesh.range(r_).far_end() : e[r_]);
    }
  }
}
