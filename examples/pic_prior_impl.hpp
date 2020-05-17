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

  struct ParticleTracing {
  public:
    static void init( const map<Properties>& properties,
                      const map<PtcArray>& particles ) {
      // initialize trace_counters to ensure unique trace serial numbers across runs
      for ( auto sp : properties ) {
        unsigned int n = 0;
        for ( const auto& ptc : particles[sp] )
          if ( ptc.is(flag::traced) )
            n = std::max<unsigned int>( n, ptc.template get<::particle::serial_number>() );

        _data()._trace_counter.insert( sp, n+1 );
      }
    }
    template < typename P >
    friend void trace( P& ptc );

  private:
    map<unsigned int> _trace_counter {}; // for assigning serial numbers to traced particles

    static ParticleTracing& _data() {
      static ParticleTracing r;
      return r;
    }

    ParticleTracing() = default;
    ParticleTracing(const ParticleTracing&);
    ParticleTracing(ParticleTracing&&) noexcept;
    ParticleTracing& operator=( const ParticleTracing& );
    ParticleTracing& operator=( ParticleTracing&& ) noexcept;
    ~ParticleTracing() = default;
  };

  template < typename P >
  void trace( P& ptc ) {
    if ( not ptc.is(flag::traced) ) {
      using sn = ::particle::serial_number;
      ptc.set(flag::traced);
      if ( ptc.template get<sn>() == 0 ) {
        ptc.set(sn(ParticleTracing::_data()._trace_counter[ptc.template get<species>()]++));
      }
    }
  }

  template < typename P >
  void untrace( P& ptc ) {
    ptc.reset(flag::traced);
  }

}

/// define convenient literal operators
constexpr pic::real_t operator"" _deg ( long double x ) noexcept {
  return static_cast<pic::real_t>(x * 3.14159265358979323846264L / 180);
}

constexpr pic::real_t operator"" _r ( long double x ) noexcept {
  return static_cast<pic::real_t>(x);
}

#include "toml++/toml.h"
using namespace std::string_view_literals;
namespace pic {
  using ConfFile_t = toml::table;

  template < typename T, typename N>
  void safe_set_from_conf( T& a, const N& node, const std::string& a_str ) {
    using U =
      std::conditional_t<std::is_same_v<T, bool>, bool,
                         std::conditional_t<std::is_floating_point_v<T>, double,
                                            std::conditional_t<std::is_integral_v<T>, int64_t, T>>>;
    auto x = node.template value<U>();
    if (!x) {
      throw std::runtime_error("ERROR : setting "+a_str+" with non-existent value in config file");
    }
    a = *x;
  }

}

#define safe_set(var,conf) safe_set_from_conf(var,conf,#var)

#endif
