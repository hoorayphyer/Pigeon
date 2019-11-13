#ifndef _IO_DATASAVER_HPP_
#define _IO_DATASAVER_HPP_

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
    std::unique_ptr<silo::file_t> master {}; // only significant on world.rank() == 0, use reference because silo library is sloppy about const-ness.
    std::string meshname;
    silo::OptList optlist;

    DataSaver( const std::optional<mpi::CartComm>& cart_opt )
      : _cart_opt(cart_opt) {}

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
        auto [c, topos] = _cart_opt->coords_topos();
        for ( int i = 0; i < topos.size(); ++i ) _block_ns_prefix += "_%03d";
        _block_ns_prefix += "/";

        _block_ns_postfix = "";
        // mpi uses row major numbering to map cartesian rank to coordinates, e.g. (0,0) -> 0, (0,1) -> 1, (1,0) -> 2, (1,1) -> 3
        std::vector<int> strides ( topos.size() + 1 ); // strides = ( DzDyDx, DzDy, Dz, 1 )
        strides.back() = 1;
        for ( int i = topos.size() - 1; i > -1; --i ) strides[i] = strides[i+1] * topos[i].dim();

        for ( int i = 0; i < topos.size(); ++i ) {
          _block_ns_postfix += delimiter + std::string("(n%" + std::to_string(strides[i]) + ")/") + std::to_string(strides[i+1]);
        }
      }
    }

    void PutMultimesh( silo::MeshType mesh_type, silo::OptList oplst ) const {
      if ( master )
        master->put_multimesh( meshname, _cart_opt->size(), _file_ns,
                               _block_ns_prefix + meshname + _block_ns_postfix,
                               mesh_type, oplst );
    }

    template < typename T, int DGrid >
    void save( std::string varname, const field::Component<T,DGrid>& fcomp ) const {
      if ( !_cart_opt ) return;
      std::vector<int> dims(DGrid);
      for ( int i = 0; i < DGrid; ++i ) dims[i] = fcomp.mesh().range(i).full_size();
      pmpio([this,varname,&fcomp,dims=std::move(dims)](auto& dbfile) {
              dbfile.put_var( varname, meshname, fcomp.data().data(), dims );
            });

      if ( master )
        master->put_multivar( varname, _cart_opt->size(), _file_ns,
                              _block_ns_prefix + varname + _block_ns_postfix,
                              optlist );
    }
  };

}

#endif
