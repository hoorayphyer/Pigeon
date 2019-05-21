#include "silopp/silo_operations.hpp"
#include <silo.h>

namespace silo {
  template <typename Type>
  constexpr int datatype(Type* = nullptr) noexcept {
    using T = std::remove_const_t<Type>;
    if constexpr ( std::is_same_v<T, char> ) return DB_CHAR;
    else if ( std::is_same_v<T, short> ) return DB_SHORT;
    else if ( std::is_same_v<T, int> ) return DB_INT;
    else if ( std::is_same_v<T, long> ) return DB_LONG;
    else if ( std::is_same_v<T, unsigned long long> ) return DB_LONG_LONG; // TODO check this
    else if ( std::is_same_v<T, float> ) return DB_FLOAT;
    else if ( std::is_same_v<T, double> ) return DB_DOUBLE;
    else return DB_NOTYPE;
  }

  template <typename T>
  constexpr int datatype(const T&) noexcept {
    return datatype((T*)0);
  }
}

namespace silo {
  namespace pmpio { struct file_t; }

  template < typename F >
  constexpr bool IsPmpio() noexcept { return false; }

  template <>
  constexpr bool IsPmpio<pmpio::file_t>() noexcept { return true; }
}

namespace silo{
  template < typename file_t >
  template < typename T >
  void Operations<file_t>::put_mesh( std::string meshname, const std::vector<std::vector<T>>& coords, MeshType mt, OptList optlist ) {
    const int ndims = coords.size();

    std::vector<const void*> raw_ptr(ndims);
    for ( int i = 0; i < ndims; ++i ) raw_ptr[i] = coords[i].data();

    std::vector<int> dims(ndims);
    for ( int i = 0; i < ndims; ++i )
      dims[i] = coords[i].size();

    auto mesh_type =
      [&mt]() {
        switch( mt ) {
        case MeshType::Curv : return DB_NONCOLLINEAR;
        default : return DB_COLLINEAR;
        }
      }();

    DBPutQuadmesh( _dbfile(), meshname.c_str(), NULL, raw_ptr.data(), dims.data(), ndims, datatype((T)0), mesh_type, optlist );
  }

  template < typename file_t >
  template < typename T >
  void Operations<file_t>::put_var( std::string varname, std::string meshname, const T* vardata, const std::vector<int>& dims, OptList optlist ) {
    DBPutQuadvar1(_dbfile(), varname.c_str(), meshname.c_str(), vardata, dims.data(), dims.size(), NULL, 0, datatype((T)0), DB_ZONECENT, optlist);
  }

  template < typename file_t >
  template < typename T >
  void Operations<file_t>::put_var( std::string varname, std::string meshname, const std::vector<const T*>& vardata, const std::vector<int>& dims, OptList optlist ) {
    int nvars = vardata.size();
    std::vector<std::string> varstrs(nvars);
    for ( int i = 0; i < nvars; ++i ) varstrs[i] = varname + std::to_string(i+1);
    std::vector<const char*> varnames(nvars);
    for ( int i = 0; i < nvars; ++i ) varnames[i] = varstrs[i].c_str();
    DBPutQuadvar(_dbfile(), varname.c_str(), meshname.c_str(), nvars, varnames.data(), vardata.data(), dims.data(), dims.size(), NULL, 0, datatype((T)0), DB_ZONECENT, optlist);
  }

  template < typename file_t >
  void Operations<file_t>:: put_multimesh( std::string multimeshname, int nblock, std::string file_ns, std::string block_ns, MeshType mt, OptList optlist ) {
    DBoptlist* raw_list = optlist;
    auto mesh_type =
      [&mt]() {
        switch( mt ) {
        case MeshType::Curv : return DB_QUAD_CURV;
        default : return DB_QUAD_RECT;
        }
      }();
    DBAddOption(raw_list, DBOPT_MB_BLOCK_TYPE, &mesh_type);
    DBAddOption(raw_list, DBOPT_MB_FILE_NS, (void*)file_ns.c_str());
    DBAddOption(raw_list, DBOPT_MB_BLOCK_NS, (void*)block_ns.c_str());
    DBPutMultimesh(_dbfile(), multimeshname.c_str(), nblock, NULL, NULL, raw_list );
  }

  template < typename file_t >
  void Operations<file_t>:: put_multivar( std::string multivarname, int nblock, std::string file_ns, std::string block_ns, OptList optlist ) {
    DBoptlist* raw_list = optlist;
    int a = DB_QUADVAR;
    DBAddOption(raw_list, DBOPT_MB_BLOCK_TYPE, &a);
    DBAddOption(raw_list, DBOPT_MB_FILE_NS, (void*)file_ns.c_str());
    DBAddOption(raw_list, DBOPT_MB_BLOCK_NS, (void*)block_ns.c_str());
    DBPutMultivar(_dbfile(), multivarname.c_str(), nblock, NULL, NULL, raw_list );
  }

  template < typename file_t >
  template < typename T >
  void Operations<file_t>::write( std::string varname, const T* vardata, const std::vector<int>& dims ) {
    DBWrite( _dbfile(), varname.c_str(), vardata, dims.data(), dims.size(), datatype((T)0) );
  }

  template < typename file_t >
  template < typename T >
  void Operations<file_t>::write( std::string varname, const std::vector<T>& vardata ) {
    // TODOL silo specifies dims as int, which may be exceeded by number of particles. Now we are relying on the fact that on one process the number of particles is smaller than int limit.
    int dims = vardata.size();
    DBWrite( _dbfile(), varname.c_str(), vardata.data(), &dims, 1, datatype((T)0) );
  }

  template < typename file_t >
  template < typename T >
  void Operations<file_t>::write( std::string varname, T vardata ) {
    int dims = 1;
    DBWrite( _dbfile(), varname.c_str(), &vardata, &dims, 1, datatype((T)0) );
  }
}

#include "silopp/silo++.hpp"
namespace silo {
  template struct Operations<file_t>;

#define INSTANTIATE_PUTTER(_TYPE_)                                      \
  template void Operations<file_t>::put_mesh( std::string, const std::vector<std::vector<_TYPE_>>&, MeshType, OptList ); \
  template void Operations<file_t>::put_var( std::string, std::string, const _TYPE_*, const std::vector<int>& dims, OptList ); \
  template void Operations<file_t>::put_var( std::string, std::string, const std::vector<const _TYPE_*>&, const std::vector<int>& dims, OptList )


#define INSTANTIATE_WRITE(_TYPE_)                                       \
  template void Operations<file_t>::write( std::string, const _TYPE_*, const std::vector<int>& dims ); \
  template void Operations<file_t>::write( std::string, const std::vector<_TYPE_>& ); \
  template void Operations<file_t>::write( std::string, _TYPE_ )


  INSTANTIATE_PUTTER(float);
  INSTANTIATE_PUTTER(double);

  INSTANTIATE_WRITE(int);
  INSTANTIATE_WRITE(double);
  INSTANTIATE_WRITE(unsigned long long);
}
