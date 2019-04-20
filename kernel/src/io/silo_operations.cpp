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
      //TODO do a downsampling outside?
    DBPutQuadvar1(_dbfile(), varname.c_str(), meshname.c_str, vardata, dims.data(), dims.size(), NULL, 0, datatype((StorageType)0), DB_NODECENT, NULL);
  }

  template < typename file_t >
  void SiloPutter<file_t>:: put_multimesh( std::string multimeshname, const std::vector<std::string>& piecenames, const OptList& optlist ) {
    std::vector<const char*> names(piecenames.size());
    for ( int i = 0; i < names.size(); ++i )
      names[i] = piecenames[i].c_str();
    DBPutMultimesh(_dbfile(), multimeshname.c_str(), names.size(), names.data(), NULL, optlist );
  }

  template < typename file_t >
  void SiloPutter<file_t>:: put_multivar( std::string multivarname, const std::vector<std::string>& piecenames, const OptList& optlist ) {
    std::vector<const char*> names(piecenames.size());
    for ( int i = 0; i < names.size(); ++i )
      names[i] = piecenames[i].c_str();
    DBPutMultivar(_dbfile(), multivarname.c_str(), names.size(), names.data(), NULL, optlist );
  }

  // auto put_multi_sth =
  //   [&pmpio=pmpio, nblocks=comm_size, isPmpio=pane.isPmpio, prefix=pane.prefix] (DBfile* dbfile, int timestep, std::string name_sth, const int type_sth, auto db_put_multisth) {
  //     // TODO move to outside
  //     // std::vector<char*> names(nblocks);

  //     // for (int i = 0; i < nblocks; i++) {
  //     //   char* name = new char[100];
  //     //   if ( isPmpio ) {
  //     //     sprintf(name, "group%d/%s%06d.d:proc%d/%s", pmpio->GroupRank(i), prefix.c_str(), timestep, i, name_sth.c_str());
  //     //   } else {
  //     //     sprintf(name, "rank%d/%s%06d.d:%s", i, prefix.c_str(), timestep, name_sth.c_str());
  //     //   }
  //     //   names.push_back(name);
  //     // }

  //     // DBoptlist* optlist = DBMakeOptlist(1);
  //     // DBAddOption(optlist, DBOPT_MB_BLOCK_TYPE, &type_sth);
  //     db_put_multisth(dbfile, name_sth.c_str(), nblocks, names.data(), NULL, optlist);
  //     // DBFreeOptlist(optlist);

  //     // for ( auto& elm : names )
  //     //   delete[] elm;
  //   };

  // put_multimesh = [=] (DBfile* dbfile, int timestep ) { return put_multi_sth(dbfile, timestep, MESHNAME, (linearMesh ? DB_QUAD_RECT : DB_QUAD_CURV), DBPutMultimesh); };


}
