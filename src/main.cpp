#include "aperture.hpp"
#include "traits.hpp"

int main() {
  Aperture< traits::real_t, traits::DGrid, traits::DPtc, traits::ptc_state_t > aperture;
  aperture.launch();
  return 0;
}
