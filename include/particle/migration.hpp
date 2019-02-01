#ifndef  _MIGRATION_HPP_
#define  _MIGRATION_HPP_

// TODO may use intercommunicator

namespace particle {
  template < typename q_t, typename borders_t >
  bool is_migrate( const q_t& q, const borders_t& bounds ) noexcept;

  template < typename ptc_array_t, typename neigh_t, typename borders_t, typename comm_t >
  struct migrate_t {
    void operator() ( ptc_array_t& buffer, const neigh_t& neighbors, const borders_t& bounds, const comm_t& comm );
  };
}

#endif
