#pragma once

namespace particle {
enum class species : int { unknown = 0, electron, positron, ion, photon };

inline constexpr int NUM_SPECIES = 5;
}  // namespace particle
inline constexpr particle::species EL = particle::species::electron;
inline constexpr particle::species PO = particle::species::positron;
inline constexpr particle::species IO = particle::species::ion;
inline constexpr particle::species PH = particle::species::photon;
