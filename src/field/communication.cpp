#include "field/communication.hpp"
#include "field/field.hpp"
#include "parallel/mpi++.hpp"
#include "apt/block.hpp"

namespace field::impl {
  template < typename Component, typename T, int D >
  void copy_to_buffer( std::vector<T>& buffer, const Component& comp,
                       const apt::Index<D>& I_b,
                       const apt::Index<D>& extent ) {
    int i = 0;
    for ( const auto& I : apt::Block(extent) )
      buffer[i++] = comp(I_b + I);
  }

}

namespace field {
  template < typename T, int DField, int DGrid >
  void sync_guard_cells( Field<T, DField, DGrid>& field, const mpi::Comm& comm ) {
    // TODOL can we use the most continous part of field for another temporary buffer? More generally, can we change the memory layout? Because guard cells in other directions are not continous anyway. Those guard cells can be used as temporary buffers.

    // Rationale: when send in direction i, include guard cells of the other
    // dimension(s). Looping order: DGrid, LeftRightness, DField

    const auto& mesh = field.mesh();
    auto I_b = mesh.origin();
    auto extent = mesh.extent();

    std::vector<T> send_buf, recv_buf;

    auto ext_size
      = []( const auto& ext ) noexcept {
          int size = 1;
          apt::foreach<0,DGrid>
            ( [&size] (auto e) { size *= e;}, ext );
          return size;
        };

    for ( int i_grid = 0; i_grid < DGrid; ++ i_grid ) {
      int ext_orig = extent[i_grid];
      extent[i_grid] = mesh.guard();
      send_buf.resize( ext_size(extent) );
      recv_buf.resize( send_buf.size() );

      int I_b_orig = I_b[i_grid];
      for ( int lr_send = 0; lr_send < 2; ++lr_send ) { // 0 is to left, 1 is to right
        I_b[i_grid] = lr_send ? ( mesh.bulk()[i_grid].dim() - mesh.guard() ) : 0;

        apt::foreach<0,DField>
          ( [&]( auto& comp )
            {
              // TODO sendrecv here.
              copy_to_buffer( send_buf, comp, I_b, extent );
              // mpi_Isend(send_buf);
              // mpi_Irecv(recv_buf);
              // mpi_wait();
              // copy_from_buffer( recv_buf, comp, I_b, extent );
            }, field );

      }
      // restroe
      I_b[i_grid] = I_b_orig;
      extent[i_grid] = ext_orig;
    }

  }
}
