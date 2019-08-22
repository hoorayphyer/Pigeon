#ifndef  _SPECIES_PREDEF_HPP_
#define  _SPECIES_PREDEF_HPP_

namespace particle {
  enum class species : char
    { unknown=0, electron, positron, ion, photon };

  inline constexpr int NUM_SPECIES = 5;
}

#endif
