#ifndef _PARTICLE_SCATTERING_HPP_
#define _PARTICLE_SCATTERING_HPP_

#include "particle/array.hpp"
#include "apt/vec.hpp"
#include <optional>
#include <memory>
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
    virtual bool operator() ( const Ptc& ) { return true; }
  };

  template < class Ptc >
  struct Channel {
    using T = typename Ptc::vec_type::element_type;
    virtual std::optional<T> operator() ( const Ptc& ptc, const apt::Vec<T,Ptc::NDim>& dp, T dt, const apt::Vec<T,Ptc::NDim>& B, util::Rng<T>& rng ) { return {}; }
  };

  template < class PtcArr >
  struct Scat {
  protected:
    using Ptc = ts::Ptc<PtcArr>;
    using T = typename Ptc::vec_type::element_type;

  private:
    std::vector<std::unique_ptr<Eligible<Ptc>>> _eligs;
    std::vector<std::unique_ptr<Channel<Ptc>>> _channels;

    virtual void impl ( std::back_insert_iterator<PtcArr> itr,
                        ts::Ptc<PtcArr>& ptc, ts::T<PtcArr> param0 ) {};

  public:
    Scat() = default;
    // NOTE deepcopy unique_ptr
    Scat( const Scat& other ) {
      _eligs.resize( other._eligs.size() );
      for ( int i = 0; i < _eligs.size(); ++i ) {
        if ( other._eligs[i] )
          _eligs[i].reset( new Eligible<Ptc>( *other._eligs[i] ) );
      }
      _channels.resize( other._channels.size() );
      for ( int i = 0; i < _channels.size(); ++i ) {
        if ( other._channels[i] )
          _channels[i].reset( new Channel<Ptc>( *other._channels[i] ) );
      }
    }


    Scat& add( const Eligible<Ptc>& elig ) {
      _eligs.emplace_back( new Eligible<Ptc>(elig) );
      return *this;
    }

    Scat& add( const Channel<Ptc>& channel ) {
      _channels.emplace_back( new Channel<Ptc>(channel) );
      return *this;
    }

    void operator() ( std::back_insert_iterator<PtcArr> itr,
                      ts::Ptc<PtcArr>& ptc, const apt::Vec<T,Ptc::NDim>& dp, T dt, const apt::Vec<T,Ptc::NDim>& B, util::Rng<ts::T<PtcArr>>& rng ) {
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
  template < bool Instant, class PtcArr >
  class RadiationFromCharges : public Scat<PtcArr> {
    virtual void impl ( std::back_insert_iterator<PtcArr> itr,
                        ts::Ptc<PtcArr>& ptc, ts::T<PtcArr> param0 ) override;
  public:
    using Scat<PtcArr>::Scat;
  };

  template < class PtcArr >
  class PhotonPairProduction : public Scat<PtcArr> {
    virtual void impl ( std::back_insert_iterator<PtcArr> itr,
                        ts::Ptc<PtcArr>& photon, ts::T<PtcArr> ) override;
  public:
    using Scat<PtcArr>::Scat;
  };
}

namespace particle {
  template < class PtcArr >
  struct ScatGen {
    void Register( species sp, const scat::Scat<PtcArr>& scat );
    void Unregister( species sp );
    scat::Scat<PtcArr>* operator() ( species sp );
  };
}

#endif
