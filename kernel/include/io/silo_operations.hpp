#ifndef  _SILO_OPERATIONS_HPP_
#define  _SILO_OPERATIONS_HPP_
#include <string>
#include <vector>

typedef struct DBfile DBfile;

namespace silo {

  struct OptList;

  template < typename file_t >
  struct SiloPutter {
  private:
    DBfile* _dbfile() noexcept { return static_cast<file_t&>(*this)._dbfile(); }

  public:
    template < typename StorageType >
    void put_mesh( std::string meshname, const std::vector<std::vector<StorageType>>& coords, const OptList& optlist );

    template < typename StorageType >
    void put_var( std::string varname, std::string meshname, const StorageType* vardata, const std::vector<int>& dims );

    void put_multimesh( std::string multimeshname, const std::vector<std::string>& piecenames, const OptList& optlist );

    void put_multivar( std::string multivarname, const std::vector<std::string>& piecenames, const OptList& optlist );
  };
}

#endif
