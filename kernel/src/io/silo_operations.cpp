#include "io/silo_operations.hpp"
#include "io/silo_optlist.hpp"
#include <silo.h>

namespace silo {
  template <typename Type>
  constexpr int datatype(Type* = nullptr) noexcept {
    using T = std::remove_const_t<Type>;
    if constexpr ( std::is_same_v<T, char> ) return DB_CHAR;
    else if ( std::is_same_v<T, short> ) return DB_SHORT;
    else if ( std::is_same_v<T, int> ) return DB_INT;
    else if ( std::is_same_v<T, long> ) return DB_LONG;
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
  template < typename StorageType >
  void SiloPutter<file_t>::put_mesh( std::string meshname, const std::vector<std::vector<StorageType>>& coords, const OptList& optlist ) {
    int ndims = coords.size();
    std::vector<int> dims(ndims);
    for (int i = 0; i < ndims; ++i )
      dims[i] = coords[i].size();

    DBPutQuadmesh(_dbfile(), meshname.c_str(), NULL,
                  coords.data(), dims.data(), ndims,
                  datatype((StorageType)0), DB_COLLINEAR, optlist);
  }

  template < typename file_t >
  template < typename StorageType >
  void SiloPutter<file_t>::put_var( std::string varname, std::string meshname, const StorageType* vardata, const std::vector<int>& dims ) {
    DBPutQuadvar1(_dbfile(), varname.c_str(), meshname.c_str, vardata, dims.data(), dims.size(), NULL, 0, datatype((StorageType)0), DB_NODECENT, NULL);
  }

  template < typename file_t >
  void SiloPutter<file_t>:: put_multimesh( std::string multimeshname, int nblock, std::string file_ns, std::string block_ns, OptList optlist ) {
    DBoptlist* raw_list = optlist;
    DBAddOption(raw_list, DBOPT_MB_BLOCK_TYPE, DB_QUAD_RECT);
    DBAddOption(raw_list, DBOPT_MB_FILE_NS, file_ns.c_str());
    DBAddOption(raw_list, DBOPT_MB_BLOCK_NS, block_ns.c_str());
    DBPutMultimesh(_dbfile(), multimeshname.c_str(), nblock, NULL, NULL, raw_list );
  }

  template < typename file_t >
  void SiloPutter<file_t>:: put_multivar( std::string multivarname, int nblock, std::string file_ns, std::string block_ns, OptList optlist ) {
    DBoptlist* raw_list = optlist;
    DBAddOption(raw_list, DBOPT_MB_BLOCK_TYPE, (void*)DB_QUADVAR);
    DBAddOption(raw_list, DBOPT_MB_FILE_NS, file_ns.c_str());
    DBAddOption(raw_list, DBOPT_MB_BLOCK_NS, block_ns.c_str());
    DBPutMultivar(_dbfile(), multivarname.c_str(), nblock, NULL, NULL, raw_list );
  }
}
