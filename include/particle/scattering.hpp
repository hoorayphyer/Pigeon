#ifndef _PARTICLE_SCATTERING_HPP_
#define _PARTICLE_SCATTERING_HPP_

#include "particle/array.hpp"
#include "apt/vec.hpp"
#include <optional>
#include "utility/rng.hpp"

namespace particle::scat::ts {
  template < class PtcArr >
  using Ptc = typename PtcArr::particle_type;

  template < class PtcArr >
  using T = typename PtcArr::value_type;
}

namespace particle::scat {
  template < class Ptc >
  struct Eligible {
    virtual bool operator() ( const Ptc& ptc ) = 0;
  };

  template < class Ptc >
  struct Channel {
    using T = typename Ptc::vec_type::element_type;
    virtual std::optional<T> operator() ( const Ptc& ptc, const apt::Vec<T,Ptc::NDim>& dp, T dt, const apt::Vec<T,Ptc::NDim>& B, util::Rng<T>& rng ) = 0;
  };

  template < class PtcArr >
  struct Scat {
  protected:
    using Ptc = ts::Ptc<PtcArr>;
    std::vector<Eligible<Ptc>*> _eligs;
    std::vector<Channel<Ptc>*> _channels;

  private:
    virtual void impl ( std::back_insert_iterator<PtcArr> itr,
                        ts::Ptc<PtcArr>& ptc, ts::T<PtcArr> param0 ) = 0;

  public:
    Scat( std::vector<Eligible<Ptc>*> eligs, std::vector<Channel<Ptc>*> channels )
      : _eligs(std::move(eligs)), _channels(std::move(channels)) {}

    void operator() ( std::back_insert_iterator<PtcArr> itr,
                      ts::Ptc<PtcArr>& ptc, util::Rng<ts::T<PtcArr>>& rng ) {
      for ( const auto elig : _eligs ) if ( !(*elig)(ptc) ) return;
      for ( const auto chnl : _channels ) {
        if ( auto param = (*chnl)(ptc, rng) ) {
          impl( itr, ptc, *param );
          return;
        }
      }
      return;
    }

  };
}

namespace particle::scat {
  template < bool Instant, class PtcArr >
  class RadiationFromCharges : private Scat<PtcArr> {
    virtual void impl ( std::back_insert_iterator<PtcArr> itr,
                        ts::Ptc<PtcArr>& ptc, ts::T<PtcArr> param0 ) override;
  public:
    using Scat<PtcArr>::Scat;
  };

  template < class PtcArr >
  class PhotonPairProduction : private Scat<PtcArr> {
    virtual void impl ( std::back_insert_iterator<PtcArr> itr,
                        ts::Ptc<PtcArr>& photon, ts::T<PtcArr> ) override;
  public:
    using Scat<PtcArr>::Scat;
  };

}

#endif
