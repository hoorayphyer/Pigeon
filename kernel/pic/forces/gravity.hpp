#ifndef _PIC_FORCE_GRAVITY_HPP_
#define _PIC_FORCE_GRAVITY_HPP_

// TODOL all the stuff under this {} are meant to be user-specified. Here the pulsar in LogSpherical is used
namespace particle::force {
  // LogSpherical
  template < typename T, template < typename > class Specs, template < typename, template < typename > class > class Ptc_t >
  void gravity( Ptc_t<T,Specs>& ptc, T dt, const apt::Vec<T,Specs<T>::Dim>& , const apt::Vec<T,Specs<T>::Dim>&, T g  ) noexcept {
    ptc.p(0) -= g * std::exp( - 2 * ptc.q(0) ) * dt;
  }
}

#endif
