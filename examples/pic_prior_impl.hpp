#ifndef _PIC_PRIOR_IMPL_HPP_
#define _PIC_PRIOR_IMPL_HPP_

#include "apt/index.hpp"
#include "apt/grid.hpp"

#include "field/field.hpp"
#include "field/haugbolle_solver/updater.hpp"

#include "particle/properties.hpp"
#include "particle/map.hpp"
#include "particle/array.hpp"
#include "particle/particle.hpp"

#include "dye/ensemble.hpp"

#include "pic.hpp"

namespace pic {
  using Index = ::apt::Index<DGrid>;
  using Grid = ::apt::Grid<real_t,DGrid>;
  using Vec3 = ::apt::Vec<real_t,3>;

  using Ensemble = ::dye::Ensemble<DGrid>;

  template < bool Const > using Component = ::field::Component<real_t,DGrid,Const>;
  using FieldAction = ::field::Action<real_t,DGrid,real_j_t>;
  using Haugbolle = ::field::Haugbolle<real_t,DGrid,real_j_t>;
  using HaugbolleBdry = ::field::HaugbolleBdry<real_t,DGrid,real_j_t>;
  template < int DField > using Field = ::field::Field<real_t,DField,DGrid>;
  using JField = ::field::Field<real_j_t,3,DGrid>;

  using ::particle::Specs;
  using ::particle::species;
  using ::particle::flag;
  using ::particle::map;
  using ::particle::Properties;
  using PtcArray = ::particle::array<real_t,Specs>;
  using PtcAction = ::particle::Action<DGrid,real_t,Specs,real_j_t>;
  using PtcUpdater = ::particle::Updater<DGrid,real_t,Specs,ShapeF,real_j_t>;
  using Particle = ::particle::Particle<real_t,Specs>;
  using Force = ::particle::Force<real_t,Specs>;
}


#endif
