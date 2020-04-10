#ifndef _MPI_COLLECTIVE_HPP_
#define _MPI_COLLECTIVE_HPP_

#include "mpipp/mpi_datatype.hpp"
#include "mpipp/mpi_request.hpp"
#include <optional>

namespace mpi {
  enum class by : char { SUM = 0, MAX, MIN, MAXLOC};

  constexpr bool IN_PLACE = true;

  template < typename Comm > // Inter can help simplify interfaces
  struct Collective_Comm {
  private:
    inline MPI_Comm _comm() const {
      return static_cast<MPI_Comm>(static_cast<const Comm&>(*this));
    }

  public:
    void barrier(const char* = "") const; // NOTE const char* is mainly for labeling barriers for users

    template < bool In_Place = false, typename T >
    std::optional<std::vector<T>>
    reduce( by op, int root, T* buffer, int count ) const;

    // TODO fix this. Returning optional on nonblocking call may not make sense
    // template < by Op, bool In_Place = false, typename T >
    // std::tuple< Request, std::optional<std::vector<T> > >
    // Ireduce( int root, T* buffer, int count ) const;

    template < typename T >
    void broadcast( int root, T* buffer, int count ) const;

    template < typename T >
    Request Ibroadcast( int root, T* buffer, int count ) const;

    template < typename T >
    std::vector<T> allgather( const T* send_buf, int send_count ) const;

    template < typename T >
    std::optional<std::vector<T>> gather( const int root, const T* send_buf, int send_count ) const;

    // scatter for intra_comm is assumed to be in_place. This means in any case, root will not send anything to itself
    template < typename T >
    void scatter( int root, T* buffer, int count ) const;

    template < typename T >
    void exscan_inplace( by op, T* send_buf, int count ) const;

    template < typename T >
    void inscan_inplace( by op, T* send_buf, int count ) const;
  };

}

#endif
