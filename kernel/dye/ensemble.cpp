#include "dye/ensemble_impl.hpp"
#include "pic.hpp"

namespace dye {
  template struct Ensemble<pic::DGrid>;

  template
  std::optional<Ensemble<pic::DGrid>> create_ensemble( const std::optional<mpi::CartComm>& cart, const std::optional<mpi::Comm>& intra );

  template
  std::optional<Ensemble<pic::DGrid>> create_ensemble( const std::optional<mpi::CartComm>& cart );
}
