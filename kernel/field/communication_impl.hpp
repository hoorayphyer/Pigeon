#include "field/communication.hpp"
#include "field/field.hpp"
#include "mpipp/mpi++.hpp"
#include "apt/block.hpp"

// TODOL can we use the most continous part of field for another temporary buffer? More generally, can we change the memory layout? Because guard cells in other directions are not continous anyway. Those guard cells can be used as temporary buffers.
namespace field {

  namespace impl {
    constexpr bool bulk_to_guard = false;
    constexpr bool guard_to_bulk = true;

    template < typename T, int DField, int DGrid, template < typename > class Policy >
    void communicate( Field<T, DField, DGrid>& field, const mpi::CartComm& comm, Policy<T> ) {
      // Looping order: DGrid, LeftRightness, DField
      const auto& mesh = field.mesh();

      auto f_Ib_send =
        [&mesh] ( int ith_dim, bool is_send_to_right ) {
          if constexpr ( Policy<T>::send_mode == bulk_to_guard ) {
              return is_send_to_right ? ( mesh.bulk_dim(ith_dim) - mesh.guard() ) : 0;
            } else {
            return is_send_to_right ? ( mesh.bulk_dim(ith_dim)  ) : - mesh.guard();
          }
        };

      auto f_Ib_recv =
        [&mesh] ( int ith_dim, bool is_send_to_right ) {
          if constexpr ( Policy<T>::send_mode == bulk_to_guard ) {
              return is_send_to_right ? - mesh.guard() : mesh.bulk_dim(ith_dim);
            } else {
            return is_send_to_right ? 0 : ( mesh.bulk_dim(ith_dim) - mesh.guard() );
          }
        };

      auto I_b = mesh.origin();
      auto extent = mesh.extent();

      std::vector<T> send_buf, recv_buf;

      for ( int i_grid = 0; i_grid < DGrid; ++ i_grid ) {
        int ext_orig = extent[i_grid];
        extent[i_grid] = mesh.guard();
        int ext_size = 1;
        apt::foreach<0,DGrid>
          ( [&ext_size] (auto e) { ext_size *= e;}, extent );

        int I_b_orig = I_b[i_grid];
        for ( int lr_send = 0; lr_send < 2; ++lr_send ) { // 0 is to left, 1 is to right

          auto [ src, dest ] = comm.shift( i_grid, lr_send ? 1 : -1 );
          // Edge case: cart dim = 1 with periodic boundaries
          if ( src && (*src == comm.rank()) ) { // NOTE check on dest is not needed
            // use I_b for I_b_send
            auto I_b_recv = I_b;
            I_b[i_grid] = f_Ib_send( i_grid, lr_send );
            I_b_recv[i_grid] = f_Ib_recv( i_grid, lr_send );
            apt::foreach<0,DField>
              ( [&]( auto comp ) { // TODOL semantics
                  for ( const auto& I : apt::Block(extent) ) {
                    T buf{}; // This is needed to accommodate merge_guard
                    Policy<T>::to_send_buf( buf, comp(I_b + I) );
                    Policy<T>::from_recv_buf( comp(I_b_recv + I), buf );
                  }
                }, field );
          } else {
            std::vector<mpi::Request> reqs(2);
            // copy all components into one buffer
            if ( dest ) {
              send_buf.resize( ext_size * DField );
              auto itr_s = send_buf.begin();
              I_b[i_grid] = f_Ib_send( i_grid, lr_send );
              apt::foreach<0,DField>
                ( [&]( auto comp ) { // TODOL semantics
                    for ( const auto& I : apt::Block(extent) )
                      Policy<T>::to_send_buf( *(itr_s++), comp(I_b + I) );
                  }, field );
              reqs[0] = comm.Isend( *dest, 924, send_buf.data(), send_buf.size() );
            }
            if ( src ) {
              recv_buf.resize( ext_size * DField );
              reqs[1] = comm.Irecv( *src, 924, recv_buf.data(), recv_buf.size() );
            }
            mpi::waitall(reqs);
            if ( src ) {
              auto itr_r = recv_buf.cbegin();
              I_b[i_grid] = f_Ib_recv( i_grid, lr_send );
              apt::foreach<0,DField>
                ( [&]( auto comp ) { // TODOL semantics
                    for ( const auto& I : apt::Block(extent) )
                      Policy<T>::from_recv_buf( comp(I_b + I), *(itr_r++) );
                  }, field );
            }
          }

        }
        // restroe
        I_b[i_grid] = I_b_orig;
        extent[i_grid] = ext_orig;
      }
    }
  }


}

namespace field {
  template < typename T >
  struct sync_policy {
    static constexpr bool send_mode = impl::bulk_to_guard;

    static constexpr void to_send_buf( T& v_buf, const T& v_data ) noexcept {
      v_buf = v_data;
    }

    static constexpr void from_recv_buf( T& v_data, const T& v_buf ) noexcept {
      v_data = v_buf;
    }
  };

  template < typename T, int DField, int DGrid >
  void sync_guard_cells_from_bulk( Field<T, DField, DGrid>& field, const mpi::CartComm& comm ) {
    impl::communicate( field, comm, sync_policy<T>{} );
  }

  template < typename T >
  struct merge_policy {
    static constexpr bool send_mode = impl::guard_to_bulk;

    static constexpr void to_send_buf( T& v_buf, T& v_data ) noexcept {
      v_buf = v_data;
      v_data = static_cast<T>(0);
    }

    static constexpr void from_recv_buf( T& v_data, const T& v_buf ) noexcept {
      v_data += v_buf;
    }
  };

  template < typename T, int DField, int DGrid >
  void merge_guard_cells_into_bulk( Field<T, DField, DGrid>& field, const mpi::CartComm& comm ) {
    impl::communicate( field, comm, merge_policy<T>{} );
  }
}
