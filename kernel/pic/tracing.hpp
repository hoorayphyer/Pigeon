#ifndef _PIC_TRACING_HPP_
#define _PIC_TRACING_HPP_

#include "particle/array.hpp"
#include "particle/map.hpp"
#include "particle/particle.hpp"
#include "particle/properties.hpp"
#include "pic/plans.hpp"

// Tracing Logic
// 1. a particle is uniquely identified by (sp, wr, sn), which we call id. wr is used in favor of ens_label because wr stays same throughout
// 2. the id is assigned only when a particle is first time traced.
// 3. sn begins with 1. sn == 0 is reversed to signal that the particle has never been traced before.
// 4. there is a flag::traced dedicated to tracing. When tracing, this flag is set, and id must be assigned accordingly. When untracing, simply reset this flag, leaving id untouched so that later another tracing doesn't assign new id to this particle
// 5. it is considered inconsistent if a particle has flag:traced set yet has sn == 0. Implementation should prevent this. In particular, setting flag::traced through state should be disabled

namespace particle {
  template <typename R, template <typename> class S>
  struct TracingManager {
  public:
    static void init( const map<Properties>& properties,
                      const map<array<R,S>>& particles ) {
      // initialize trace_counters to ensure unique trace serial numbers across runs
      for ( auto sp : properties ) {
        unsigned int n = 0;
        for ( const auto& ptc : particles[sp] )
          if ( ptc.is(flag::traced) )
            n = std::max<unsigned int>( n, ptc.template get<pid>() );

        _data()._trace_counter.insert( sp, n+1 );
      }
    }

    template < typename P >
    static void trace( P& ptc ) {
      if (not ptc.is(flag::traced)) {
        ptc.set(flag::traced);
        if (ptc.template get<pid>() == 0) {
          ptc.set(sn(_data()._trace_counter[ptc.template get<species>()]++));
        }
      }
    }

    template <typename P>
    static void untrace(P &ptc) { ptc.reset(flag::traced); }

  private:
    unsigned int _wr {}; // world rank
    map<unsigned int> _trace_counter {}; // for assigning serial numbers to traced particles

    static TracingManager& _data() {
      static TracingManager r;
      return r;
    }

    TracingManager() = default;
    TracingManager(const TracingManager&);
    TracingManager(TracingManager&&) noexcept;
    TracingManager& operator=( const TracingManager& );
    TracingManager& operator=( TracingManager&& ) noexcept;
    ~TracingManager() = default;
  };

}

namespace particle {
  template <int DGrid, typename R, template <typename> class S, typename RJ>
  struct Tracer : public Action<DGrid,R,S,RJ> {
  private:
    std::vector<species> _sps;

    bool _is_check_within_range = true;

    pic::Plan _plan{};

    using FMark_t = void (*)(typename array<R,S>::particle_type &ptc, util::Rng<R> &rng);
    FMark_t _marker = nullptr;

  public:
    Tracer* Clone() const override {return new auto(*this);}

    auto& set_marker( FMark_t f ) noexcept { _marker = f; return *this;}
    auto& set_species(const std::vector<species>& sps) noexcept {
      _sps = sps; return *this;
    }
    auto& set_is_check_within_range( bool a ) noexcept {_is_check_within_range=a; return *this;}
    auto& set_plan( const pic::Plan& p ) noexcept { _plan = p; return *this; }

    static bool is_within_bounds(const typename array<R,S>::particle_type::vec_type &q,
                                 const apt::array<apt::array<R, 2>, DGrid> &bds) {
      for (int i = 0; i < DGrid; ++i) {
        if (q[i] < bds[i][0] or q[i] >= bds[i][1])
          return false;
      }
      return true;
    }

    void operator() ( map<array<R,S>>& particles, field::Field<RJ,3,DGrid>& ,
                      std::vector<Particle<R,S>>* ,
                      const map<Properties>& ,
                      const field::Field<R,3,DGrid>&, const field::Field<R,3,DGrid>&,
                      const apt::Grid<R,DGrid>& grid, const dye::Ensemble<DGrid>* ,
                      R , int timestep, util::Rng<R>& rng) override {
      if ( !_plan.is_do(timestep) or apt::range::is_empty(*this) or !_marker ) return;

      apt::array< apt::array<R,2>, DGrid > bds;
      for ( int i = 0; i < DGrid; ++i ) {
        bds[i][0] = grid[i].absc( apt::range::begin(*this,i) );
        bds[i][1] = grid[i].absc( apt::range::end(*this,i) );
      }

      for ( auto sp : _sps ) {
        for ( auto ptc : particles[sp] ) { // TODOL semantics
          if ( !ptc.is(flag::exist)
               or (_is_check_within_range and !is_within_bounds(ptc.q(), bds) )
               ) continue;
          _marker(ptc,rng);
        }
      }
    }
  };
}

#endif
