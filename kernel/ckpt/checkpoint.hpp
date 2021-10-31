#pragma once

#include <optional>
#include <string>

#include "dye/ensemble.hpp"

namespace field {
template <typename, int, int>
struct Field;
}

namespace particle {
template <typename>
struct map;
template <typename, template <typename> class>
struct array;
struct Properties;
}  // namespace particle

// TODO other things to save, like initial conditions, shape functions? No, use
// gen.hpp and build a target called resume
namespace ckpt {
template <int DGrid, typename R, template <typename> class S>
std::string save_checkpoint(
    std::string prefix, const int num_parts,
    const std::optional<dye::Ensemble<DGrid>>& ens_opt, int timestep,
    const field::Field<R, 3, DGrid>& E, const field::Field<R, 3, DGrid>& B,
    const particle::map<particle::array<R, S>>& particles,
    const particle::map<particle::Properties>& properties);

template <int DGrid, typename R, template <typename> class S>
int load_checkpoint(std::string dir,
                    std::optional<dye::Ensemble<DGrid>>& ens_opt,
                    const std::optional<mpi::CartComm>& cart_opt,
                    field::Field<R, 3, DGrid>& E, field::Field<R, 3, DGrid>& B,
                    particle::map<particle::array<R, S>>& particles,
                    const particle::map<particle::Properties>& properties,
                    int target_load = 0);
}  // namespace ckpt
