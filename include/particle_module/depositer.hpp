#ifndef _DEPOSITER_HPP_
#define _DEPOSITER_HPP_

// TODO check missing 4pi's maybe
namespace particle {
  // NOTE Trl here may use long double
  template < sf::shape S, typename T_WJ, typename Tvt, std::size_t DPtc, std::size_t DField,
             std::size_t DGrid, typename Trl = vec::remove_cvref_t<Tvt> >
  void depositWJ( Field<T_WJ,DField>& WJ, const Particle<Tvt,DPtc>& ptc, const Vec<Trl,DPtc>& dq, const Grid<DGrid, Trl>& grid );
}

#endif
