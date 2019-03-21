#ifndef _IO_EXPORTEE_HPP_
#define _IO_EXPORTEE_HPP_

#include "apt/index.hpp"
#include "particle/load_type.hpp"

namespace field { template < typename, int > struct Mesh; }

namespace io {
  template < typename T, int DGrid >
  struct FieldBasedExportee {
    virtual T operator() ( const apt::Index<DGrid>& I_export_mesh_bulk,
                           const field::Mesh<T,DGrid>& export_mesh ) = 0;
  };
}

namespace io {
  template < typename T, int DGrid >
  struct ParticleBasedExportee {
    using particle::load_t;
    virtual load_t begin() const = 0;
    virtual load_t end() const = 0;
    load_t& next( load_t& index ) const = 0; // used for skipping empty particles for example, or iterate from back
    virtual T val() ( load_t index ) const = 0;
    virtual apt::array<T, DGrid> loc() ( load_t index ) const = 0;
  };
}

#endif
