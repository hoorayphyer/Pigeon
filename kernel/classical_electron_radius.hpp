#ifndef _CLASSICAL_ELECTRON_RADIUS_
#define _CLASSICAL_ELECTRON_RADIUS_

namespace pic {
  constexpr real_t classic_electron_radius () noexcept {
    real_t res = wdt_pic * wdt_pic / ( 4 * std::acos(-1.0l) * dt * dt);
    apt::foreach<0,DGrid>( [&res](const auto& g) { res *= g.delta(); }, supergrid );
    return res;
  }
}

#endif
