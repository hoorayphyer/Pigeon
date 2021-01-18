#include "field/log_spherical_solver/updater.hpp"
#include <optional>
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
        _fr[i] = std::exp(-2*grid[r_].absc(b_r+i));

      auto recip_sin = [thr = 0.5 * 0.5 * grid[th_].delta()]( auto th ) {
        auto x = std::sin(th);
        x = ( std::abs(x) < thr ) ? 0.0 : 1.0 / x;
        return x;
      };

      for (int j = 0; j < e_th - b_th; ++j) {
        auto th = grid[th_].absc(b_th + j);
        _sin_th[2 * j] = recip_sin(th);
        _sin_th[2 * j + 1] = 1.0 / std::sin(th + 0.5 * grid[th_].delta());
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
  struct Diff_r {
    const _helper<R,D>& h;
    const R delta;

    void operator()(Component<R,D,false>&& output, const Field<R, 3, D> &F, const bool is_E_type, const int comp, R prefactor, const int (&b_in)[D], const int (&e_in)[D], const bool lo, const bool hi) { // TODOL sematics
      prefactor /= delta;

      const auto& m = F.mesh();
      const auto s = m.stride(r_);

      if (is_E_type and lo) { // conductor surface
        const auto g1 = h.estag() * h.estag();
        auto f = prefactor * h.estag();
        for (int j = b_in[th_]; j < e_in[th_]; ++j) {
          int li = m.linear_index(th_, j) + b_in[r_] * s;
          auto x = F[comp][li + s] * g1 - F[comp][li];
          output[li] +=
              f * h.fr(b_in[r_]) *
              (x + x + g1 * (F[comp][li + s] - F[comp][li + (s << 1)] * g1));
        }
      }

      { // bulk
        R f = is_E_type ? prefactor * h.estag() : prefactor;
        R g = 1.0 / (h.estag() * h.estag());
        for (int j = b_in[th_]; j < e_in[th_]; ++j) {
          int li_j = m.linear_index(th_, j);
          for (int i = b_in[r_] + is_E_type; i < e_in[r_] - !is_E_type; ++i) {
            int li = li_j + i * s;
            output[li] += f * h.fr(i) *
                          (F[comp][li + (!is_E_type) * s] -
                           F[comp][li - (is_E_type)*s] * g);
          }
        }
      }

      // at r outer boundary, the r derivative is taken to vanish.
    }
  };

  template < typename R, int D >
  struct Diff_th {
    const _helper<R, D> &h;
    const R delta;

    void operator()(Component<R,D,false>&& output, const Field<R, 3, D> &F, const bool is_E_type, const int comp, R prefactor, const int (&b_in)[D], const int (&e_in)[D], const bool lo, const bool hi) {// TODOL sematics
      prefactor /= delta;
      const auto& m = F.mesh();
      const auto s = m.stride(th_);

      if ( comp == r_ ) {
        const R f = is_E_type ? prefactor : prefactor / h.estag();
        for (int j = b_in[th_]+is_E_type; j < e_in[th_] - !is_E_type; ++j) {
          int li_j = m.linear_index(th_, j);
          for (int i = b_in[r_]; i < e_in[r_]; ++i) {
            int li = li_j + i * m.stride(r_);
            output[li] += f * h.fr(i) * (F[comp][li + (!is_E_type)*s] - F[comp][li-(is_E_type)*s]);
          }
        }
      }
      else if (!is_E_type) {
        for (int j = b_in[th_]; j < e_in[th_] - 1; ++j) {
          auto f = prefactor / h.fsin(j,true);
          int li_j = m.linear_index(th_, j);
          for (int i = b_in[r_]; i < e_in[r_]; ++i) {
            int li = li_j + i * m.stride(r_);
            output[li] += f * h.fr(i) * (F[comp][li + s] * h.fsin(j+1) - F[comp][li] * h.fsin(j));
          }
        }
      }
      else {
        const R f = prefactor / h.estag();
        auto diff_bdry
          = [&output, e = e_in[th_], &h=this->h, &F, &comp, b_r=b_in[r_], e_r=e_in[r_], &f]( const bool lo ) {
            const auto& m = F.mesh();

            const R ff = lo ? 4.0 * f : -4.0 * f;
            const int y = lo ? 0 : -1;

            const int li_j = m.linear_index(th_, lo ? 0 : e);
            for (int i = b_r; i < e_r; ++i) {
              int li = li_j + i * m.stride(r_);
              // NOTE order O(dth^2)
              output[li] += ff * h.fr(i) * F[comp][li + y * m.stride(th_)];
            }
          };

        if (lo) diff_bdry(true);

        for (int j = b_in[th_]+1; j < e_in[th_]; ++j) {
          int li_j = m.linear_index(th_, j);
          R fsr = f * h.fsin(j, true) / h.fsin(j);
          R sr = h.fsin(j - 1, true) / h.fsin(j, true);
          for (int i = b_in[r_]; i < e_in[r_]; ++i) {
            int li = li_j + i * m.stride(r_);
            output[li] += fsr * h.fr(i) * (F[comp][li] - F[comp][li - m.stride(th_)] * sr);
          }
        }

        if (hi) diff_bdry(false);
      }

    }

  };
}

namespace field {
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

    _helper<R, DGrid> h{grid, b[r_], mesh.range(r_).far_end(), b[th_], e[th_]}; // NOTE e[r_] is not used here

    Diff_r<R,DGrid> diff_r {h, grid[r_].delta()};
    Diff_th<R, DGrid> diff_th{h, grid[th_].delta()};

    auto mul_absorb =
      []( auto& output, const auto& input, auto factor ){
        // assume input and output have same mesh
        const auto& m = input.mesh();
        for ( int c = 0; c < 3; ++c ) {
          for ( std::size_t i = 0; i < input[c].data().size(); ++i )
            output[c].data()[i] += factor * input[c].data()[i];
        }
      };

    auto curl_E
      = [&diff_r, &diff_th,lo=bool(lo_axis),hi=bool(hi_axis),surf=bool(surf), outer=bool(outer), &b, &e ](auto &F, const auto &E, const auto factor, bool is_continuous = false) mutable {
      e[th_] -= hi;
      diff_th(F[r_], E, true, phi_, factor, b, e, lo, hi);
      if ( surf and is_continuous ) b[r_]--;
      diff_r(F[th_], E, true, phi_, -factor, b, e, surf and not is_continuous, outer);
      e[th_] += hi;
      diff_r(F[phi_], E, true, th_, factor, b, e, surf and not is_continuous, outer); // NOTE d_r E_th is NOT enforced to be zero on higher axis
      if (surf and is_continuous) b[r_]++;
      e[th_] -= hi;
      diff_th(F[phi_], E, true, r_, -factor, b, e, lo, hi);
      e[th_] += hi;

      if (!surf) b[r_]++;
      if (!lo) b[th_]++;
    };

    auto curl_B
      = [&diff_r, &diff_th,lo=bool(lo_axis),hi=bool(hi_axis),surf=bool(surf), outer=bool(outer), &b, &e](auto &F, const auto &B, const auto factor) mutable {
      diff_th(F[r_], B, false, phi_, factor, b, e, lo, hi);
      diff_r(F[th_], B, false, phi_, -factor, b, e, surf, outer); // NOTE d_r B_phi is NOT enforced to be zero on higher axis
      e[th_] -= hi;
      diff_r(F[phi_], B, false, th_, factor, b, e, surf, outer);
      e[th_] += hi;
      diff_th(F[phi_], B, false, r_, -factor, b, e, lo, hi);

      if (!outer) e[r_]--;
      if (!hi) e[th_]--;
    };

    {
      for (int j = b[th_]; j < e[th_] - bool(hi_axis); ++j) {
        int li_j = mesh.linear_index(th_, j);
        R fsr =  -_fourpi * dt * h.fsin(j,true);
        for (int i = b[r_]; i < e[r_]; ++i) {
          int li = li_j + i * mesh.stride(r_);
          J[0][li] *= fsr * h.fr(i);
        }
      }
      for (int j = b[th_]; j < e[th_]; ++j) {
        int li_j = mesh.linear_index(th_, j);
        R fsr = -_fourpi * dt * h.fsin(j) / ( h.estag() * h.estag());
        for (int i = b[r_]; i < e[r_] - bool(outer); ++i) {
          int li = li_j + i * mesh.stride(r_);
          J[1][li] *= fsr * h.fr(i);
        }
      }
      for (int j = b[th_]; j < e[th_] - bool(hi_axis); ++j) {
        int li_j = mesh.linear_index(th_, j);
        R fsr =  -_fourpi * dt / (h.estag() * h.estag());
        for (int i = b[r_]; i < e[r_]-bool(outer); ++i) {
          int li = li_j + i * mesh.stride(r_);
          J[2][li] *= fsr * h.fr(i);
        }
      }

      for (auto &x : h.fr()) x = std::sqrt(x);
      for (auto &x : h.fsin()) {
        x = (std::abs(x) < 0.5 * 0.5 * grid[th_].delta()) ? 0.0 : 1.0 / x;
      }
    }

    { // B version
      curl_B(J, B, (1.0 - alpha)*dt);
      curl_E(B, J, -alpha*dt);
      curl_E(B, E, -dt,  true);
      mul_absorb(E, J, 1.0);
      auto C = B;

      for (int i = std::max(_op_inv_precision, 1); i != 1; --i) {
        J.reset();
        curl_B(J, C, 1.0);
        C.reset();
        curl_E(C, J, -(alpha * alpha * dt * dt));
        mul_absorb(B, C, 1.0);
      }
      J.reset();
      curl_B(J, C, 1.0);
      curl_E(B, J, -(alpha * alpha * dt * dt));

      curl_B(E, B, (alpha * dt));
    }
    // { // E version
    //   Field<R, 3, DGrid> C(mesh);
    //   curl_E(C, E, 1.0, true);

    //   mul_absorb(B, C, -alpha * (1-alpha) * dt);
    //   curl_B(J, B, dt);
    //   mul_absorb(B, C, -(1 - alpha) * (1 - alpha) * dt);

    //   mul_absorb(E, J, 1.0);
    //   for (int i = std::max(_op_inv_precision, 1); i != 1; --i) {
    //     curl_E(C, J, 1.0);
    //     J.reset();
    //     curl_B(J, C, -(alpha * alpha * dt * dt));
    //     C.reset();
    //     mul_absorb(E, J, 1.0);
    //   }
    //   curl_E(C, J, 1.0);
    //   curl_B(E, C, -(alpha * alpha * dt * dt));

    //   curl_E(B, E, -(alpha * dt));
    // }
    // NOTE sync E B is done outside
  }
}
