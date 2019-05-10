#ifndef _PARTICLE_SCATTERING_HPP_
#define _PARTICLE_SCATTERING_HPP_

#include "particle/array.hpp"
#include "apt/vec.hpp"
#include <optional>
#include <memory>
#include "random/rng.hpp"

namespace particle::scat {
  template < typename T, template < typename > class PtcSpecs >
  using Ptc_t = typename array<T,PtcSpecs>::particle_type;

  template < typename T, template < typename > class PtcSpecs >
  struct Eligible {
    virtual bool operator() ( const Ptc_t<T,PtcSpecs>& ) { return true; }
  };

  template < typename T, template < typename > class PtcSpecs >
  struct Channel {
    virtual std::optional<T> operator() ( const Ptc_t<T,PtcSpecs>& ptc, const apt::Vec<T,PtcSpecs<T>::Dim>& dp, T dt, const apt::Vec<T,PtcSpecs<T>::Dim>& B, util::Rng<T>& rng ) { return {}; }
  };

  template < typename T, template < typename > class PtcSpecs >
  struct Scat {
  protected:
    using Ptc = Ptc_t<T,PtcSpecs>;

  private:
    std::vector<std::unique_ptr<Eligible<T,PtcSpecs>>> _eligs;
    std::vector<std::unique_ptr<Channel<T,PtcSpecs>>> _channels;

    virtual void impl ( std::back_insert_iterator<array<T,PtcSpecs>> itr,
                        Ptc& ptc, T param0 ) {};

  public:
    Scat() = default;
    // NOTE deepcopy unique_ptr
    Scat( const Scat& other ) {
      _eligs.resize( other._eligs.size() );
      for ( int i = 0; i < _eligs.size(); ++i ) {
        if ( other._eligs[i] )
          _eligs[i].reset( new Eligible( *other._eligs[i] ) );
      }
      _channels.resize( other._channels.size() );
      for ( int i = 0; i < _channels.size(); ++i ) {
        if ( other._channels[i] )
          _channels[i].reset( new Channel( *other._channels[i] ) );
      }
    }


    Scat& add( const Eligible<T,PtcSpecs>& elig ) {
      _eligs.emplace_back( new Eligible(elig) );
      return *this;
    }

    Scat& add( const Channel<T,PtcSpecs>& channel ) {
      _channels.emplace_back( new Channel(channel) );
      return *this;
    }

    void operator() ( std::back_insert_iterator<array<T,PtcSpecs>> itr,
                      Ptc& ptc, const apt::Vec<T,PtcSpecs<T>::Dim>& dp, T dt, const apt::Vec<T,PtcSpecs<T>::Dim>& B, util::Rng<T>& rng ) {
      for ( const auto& elig : _eligs ) if ( !(*elig)(ptc) ) return;
      for ( const auto& chnl : _channels ) {
        if ( auto param = (*chnl)(ptc, dp, dt, B, rng) ) {
          impl( itr, ptc, *param );
          return;
        }
      }
      return;
    }

  };
}

namespace particle::scat {
  template < bool Instant, typename T, template < typename > class PtcSpecs >
  class RadiationFromCharges : public Scat<T,PtcSpecs> {
    virtual void impl ( std::back_insert_iterator<array<T,PtcSpecs>> itr,
                        Ptc_t<T,PtcSpecs>& ptc, T param0 ) override;
  public:
    using Scat<T,PtcSpecs>::Scat;
  };

  template < typename T, template < typename > class PtcSpecs >
  class PhotonPairProduction : public Scat<T,PtcSpecs> {
    virtual void impl ( std::back_insert_iterator<array<T,PtcSpecs>> itr,
                        Ptc_t<T,PtcSpecs>& photon, T ) override;
  public:
    using Scat<T,PtcSpecs>::Scat;
  };
}

namespace particle {
  template < typename T, template < typename > class PtcSpecs >
  struct ScatGen {
    void Register( species sp, const scat::Scat<T,PtcSpecs>& scat );
    void Unregister( species sp );
    scat::Scat<T,PtcSpecs>* operator() ( species sp );
  };
}

#endif
