#ifndef _IO_AUX_HPP_
#define _IO_AUX_HPP_

namespace io {
  constexpr int POW( int B, int E ) {
    if ( E == 0 ) return 1;
    else return B * POW(B,E-1);
  };

  // TODO
  // template < typename T, int DGrid >
  // constexpr auto q_from_cell( const apt::Index<DGrid>& I, const apt::array<field::offset_t, DGrid>& ofs ) noexcept {
  //   apt::array<T,DGrid> res;
  //   for ( int i = 0; i < DGrid; ++i ) res[i] = I[i] + static_cast<T>(ofs[i]);
  //   return res;
  // }

  template < typename T, int DGrid >
  class Downsampler {
  private:
    int _ratio = 1;
    field::Field<T,1,DGrid> _fds; // downsampled field
    const std::optional<mpi::CartComm>& _cart_opt;

  public:
    Downsampler( int ratio, apt::Index<DGrid> full_bulk_dims, int guard_of_sampled_field, const std::optional<mpi::CartComm>& cart_opt )
      : _ratio(ratio), _cart_opt(cart_opt) {
      for ( auto& x: full_bulk_dims ) x /= _ratio;
      _fds = {{ full_bulk_dims, guard_of_sampled_field }};
      // _fds is defined to have MIDWAY in all directions
      for ( int i = 0; i < DGrid; ++i ) _fds.set_offset(0,i,MIDWAY);
    }

    // TODO This doesn't interpolate to zone center yet. Also may need intermediate sync_guards when interpolating to center. NOTE downsampling interpolation has its own shape function
    template < typename U >
    const field::Field<T,1,DGrid>& operator() ( const field::Component<U,DGrid>& fcomp ) {
      apt::Index<DGrid> subext;
      for ( int i = 0; i < DGrid; ++i ) subext[i] = _ratio;

      for ( const auto& Ids : apt::Block( _fds.mesh().bulk_dims() ) ) {
        U f = 0.0;
        for ( const auto& Isub : apt::Block(subext) ) {
          apt::Index<DGrid> I;
          for ( int i = 0; i < DGrid; ++i )
            I[i] = _ratio * Ids[i] + Isub[i];
          f += fcomp( I ); // TODO use interpolation
        }
        _fds[0](Ids) = f / POW(_ratio, DGrid);
      }
      if ( _cart_opt )
        field::copy_sync_guard_cells( _fds, *_cart_opt );
      return _fds;
    }

  };

}

namespace io {
  class DataSaver {
  private:
    const std::optional<mpi::CartComm>& _cart_opt;

    // following are significant only on cart.rank() == 0
    std::string _file_ns;
    std::string _block_ns_prefix;
    std::string _block_ns_postfix;

  public:
    silo::Pmpio pmpio; // only significant on carts
    silo::file_t master; // only significant on world.rank() == 0
    std::string meshname;
    silo::OptList optlist;

    DataSaver( const std::optional<mpi::CartComm>& cart_opt ) : _cart_opt(cart_opt) {}

    void set_namescheme(std::string str_ts, int num_files) {
      if ( !_cart_opt ) return;

      // set up file_ns and block_ns. n below is thought of as the cartesian rank. Only significant on world.rank() == 0
      constexpr char delimiter = '|';
      // NOTE file namescheme uses relative path so that when the directory is moved, the data is still viewable
      {
        _file_ns = delimiter + std::string("data/timestep") + str_ts + "/set%d.silo" + delimiter + "n%" + std::to_string(num_files);
      }

      {
        _block_ns_prefix = delimiter + std::string("cart");
        auto [c, dims, p] = _cart_opt->coords_dims_periodic();
        for ( int i = 0; i < dims.size(); ++i ) _block_ns_prefix += "_%03d";
        _block_ns_prefix += "/";

        _block_ns_postfix = "";
        // mpi uses row major numbering to map cartesian rank to coordinates, e.g. (0,0) -> 0, (0,1) -> 1, (1,0) -> 2, (1,1) -> 3
        std::vector<int> strides ( dims.size() + 1 ); // strides = ( DzDyDx, DzDy, Dz, 1 )
        strides.back() = 1;
        for ( int i = dims.size() - 1; i > -1; --i ) strides[i] = strides[i+1] * dims[i];

        for ( int i = 0; i < dims.size(); ++i ) {
          _block_ns_postfix += delimiter + std::string("(n%" + std::to_string(strides[i]) + ")/") + std::to_string(strides[i+1]);
        }
      }
    }

    void PutMultimesh( auto mesh_type, silo::OptList oplst ) {
      if ( _cart_opt && _cart_opt->rank() == 0 )
        master.put_multimesh( meshname, _cart_opt->size(), _file_ns,
                              _block_ns_prefix + meshname + _block_ns_postfix,
                              mesh_type, oplst );
    }

    void PutMultivar(std::string varname ) {
      if ( _cart_opt && _cart_opt->rank() == 0 )
        master.put_multivar( varname, _cart_opt->size(), _file_ns,
                             _block_ns_prefix + varname + _block_ns_postfix,
                             optlist );
    }

    template < typename T, int DGrid >
    void save( std::string varname, const field::Field<T,1,DGrid>& field ) {
      if ( !_cart_opt ) return;
      std::vector<int> dims(DGrid);
      for ( int i = 0; i < DGrid; ++i ) dims[i] = field.mesh().extent()[i];
      pmpio([this,varname,&field,dims=std::move(dims)](auto& dbfile) {
              dbfile.put_var( varname, meshname, field[0].data().data(), dims );
            });
      PutMultivar(varname);
    }
  };
}

#endif
