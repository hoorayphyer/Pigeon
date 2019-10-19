#include "field/sync.hpp"
#include "field/field.hpp"
#include "mpipp/mpi++.hpp"
#include "apt/block.hpp"
#include <cassert>

// TODOL can we use the most continous part of field for another temporary buffer? More generally, can we change the memory layout? Because guard cells in other directions are not continous anyway. Those guard cells can be used as temporary buffers.
namespace field {

  namespace impl {
    template < typename T, int DField, int DGrid, template < typename > class Policy >
    void sync( Field<T, DField, DGrid>& field, const mpi::CartComm& comm,
               const apt::array<apt::Range,DGrid>& domain_range,
               const apt::array<apt::Range,DGrid>& sub_range, Policy<T> ) {
      // Looping order ( from outer to inner ): DGrid, LeftRightness, DField
      assert(apt::range::is_margin_uniform(domain_range)); // POLEDANCE this implementation requires all margins to be equal
      // check sub_range \in field.range \in domain_range
      {
        const auto& range = field.mesh().range();
        bool ordered = true;
        for ( int i = 0; i < DGrid; ++i ) {
          bool x = (sub_range[i].far_begin() >= range[i].far_begin() and range[i].far_begin() >= domain_range[i].far_begin() );
          ordered = ( ordered and x );
          x = (sub_range[i].far_end() <= range[i].far_end() and range[i].far_end() >= domain_range[i].far_end() );
          ordered = ( ordered and x );
        }
        assert(ordered); // POLEDANCE print error message
      }

      const int guard = domain_range[0].margin()[LFT];

      std::vector<T> send_buf, recv_buf;

      for ( int i_grid = 0; i_grid < DGrid; ++ i_grid ) {
        { // determine if communication is needed in this direction. Each dimension has two sides, and each side must be fully contained ( i.e. full_size() ) in domain_range bulk to avoid communication, unless that side corresponds to a true boundary.
          auto neigh = comm.shift( i_grid, 1 );
          bool in_bulk = true;
          {
            bool is_contained = (sub_range[i_grid].far_begin() >= domain_range[i_grid].begin() );
            in_bulk = (in_bulk and ( is_contained or neigh[LFT] ) );
            is_contained = (sub_range[i_grid].far_end() <= domain_range[i_grid].end() );
            in_bulk = (in_bulk and ( is_contained or neigh[RGT] ) );
          }
          if ( in_bulk ) continue;
        }

        auto I_b = apt::range::far_begin(sub_range);
        auto extent = apt::range::full_size(sub_range);

        extent[i_grid] = Policy<T>::extent(guard);
        int ext_size = 1;
        for ( int i = 0; i < DGrid; ++i ) ext_size *= extent[i];

        for ( int lr_send = 0; lr_send < 2; ++lr_send ) { // 0 is to left, 1 is to right

          auto [ src, dest ] = comm.shift( i_grid, lr_send ? 1 : -1 );
          // Edge case: cart dim = 1 with periodic boundaries POLEDANCE see front
          if ( src && (*src == comm.rank()) ) { // NOTE check on dest is not needed
            // use I_b for I_b_send
            auto I_b_recv = I_b;
            I_b[i_grid] = Policy<T>::sendIb( sub_range[i_grid].size(), guard, lr_send );
            I_b_recv[i_grid] = Policy<T>::recvIb( sub_range[i_grid].size(), guard, lr_send );
            for ( int C = 0; C < DField; ++ C ) {
              for ( const auto& I : apt::Block({}, extent) ) {
                T buf{}; // This is needed to accommodate merge_guard
                Policy<T>::to_send_buf( buf, field[C](I_b + I) );
                Policy<T>::from_recv_buf( field[C](I_b_recv + I), buf, lr_send );
              }
            }
          } else {
            std::vector<mpi::Request> reqs(2);
            // copy all components into one buffer
            if ( dest ) {
              send_buf.resize( ext_size * DField );
              auto itr_s = send_buf.begin();
              I_b[i_grid] = Policy<T>::sendIb( sub_range[i_grid].size(), guard, lr_send );
              for ( int C = 0; C < DField; ++C ) {
                for ( const auto& I : apt::Block(I_b, I_b + extent) )
                  Policy<T>::to_send_buf( *(itr_s++), field[C](I) );
              }
              reqs[0] = comm.Isend( *dest, 924, send_buf.data(), send_buf.size() );
            }
            if ( src ) {
              recv_buf.resize( ext_size * DField );
              reqs[1] = comm.Irecv( *src, 924, recv_buf.data(), recv_buf.size() );
            }
            mpi::waitall(reqs);
            if ( src ) {
              auto itr_r = recv_buf.cbegin();
              I_b[i_grid] = Policy<T>::recvIb( sub_range[i_grid].size(), guard, lr_send );
              for ( int C = 0; C < DField; ++C ) {
                for ( const auto& I : apt::Block(I_b, I_b + extent) )
                  Policy<T>::from_recv_buf( field[C](I), *(itr_r++), lr_send );
              }
            }
          }
        }
      }
    }
  }


}

namespace field {
  template < typename T >
  struct copy_policy {
    static constexpr bool is_copy_sync = true;

    static constexpr void to_send_buf( T& v_buf, const T& v_data ) noexcept {
      v_buf = v_data;
    }

    static constexpr void from_recv_buf( T& v_data, const T& v_buf, bool ) noexcept {
      v_data = v_buf;
    }

    static constexpr int sendIb( int bulk_dim, int guard, bool is_send_to_right ) noexcept {
      return is_send_to_right ? ( bulk_dim - guard ) : 0;
    }

    static constexpr int recvIb( int bulk_dim, int guard, bool is_send_to_right ) noexcept {
      return is_send_to_right ? - guard : bulk_dim;
    }

    static constexpr int extent( int guard ) noexcept { return guard; }
  };

  template < typename T, int DField, int DGrid >
  void copy_sync_guard_cells( Field<T, DField, DGrid>& field, const mpi::CartComm& comm,
                              const apt::array<apt::Range,DGrid>& domain_range,
                              const apt::array<apt::Range,DGrid>& sub_range ) {
    impl::sync( field, comm, domain_range, sub_range, copy_policy<T>{} );
  }

  template < typename T >
  struct merge_policy {
    static constexpr bool is_copy_sync = false;

    static constexpr void to_send_buf( T& v_buf, T& v_data ) noexcept {
      v_buf = v_data;
    }

    static constexpr void from_recv_buf( T& v_data, const T& v_buf, bool is_send_to_right ) noexcept {
      if ( is_send_to_right )
        v_data = v_buf;
      else
        v_data += v_buf;
    }

    static constexpr int sendIb( int bulk_dim, int guard, bool is_send_to_right ) noexcept {
      return is_send_to_right ? bulk_dim - guard : - guard;
    }

    static constexpr int recvIb( int bulk_dim, int guard, bool is_send_to_right ) noexcept {
      return is_send_to_right ? -guard : ( bulk_dim - guard );
    }

    static constexpr int extent( int guard ) noexcept { return 2 * guard; }
  };

  template < typename T, int DField, int DGrid >
  void merge_sync_guard_cells( Field<T, DField, DGrid>& field, const mpi::CartComm& comm,
                               const apt::array<apt::Range,DGrid>& domain_range,
                               const apt::array<apt::Range,DGrid>& sub_range ) {
    impl::sync( field, comm, domain_range, sub_range, merge_policy<T>{} );
  }
}
