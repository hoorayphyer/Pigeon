#ifndef  _MIGRATION_HPP_
#define  _MIGRATION_HPP_

// TODO may use intercommunicator

namespace particle {
  template < typename T_q, typename T_bounds >
  bool is_migrate( const T_q& q, const T_bounds& bounds ) noexcept;

  template < typename PtcVector, typename T_neigh, typename T_bounds, typename Comm >
  void migrate( PtcVector& buffer, const T_neigh& neighbors, const T_bounds& bounds, const Comm& comm );
}

#endif
