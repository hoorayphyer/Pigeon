#ifndef  _SILO_OPTLIST_HPP_
#define  _SILO_OPTLIST_HPP_

#include <variant>
#include <unordered_map>
#include <memory>

typedef struct DBoptlist_ DBoptlist;

namespace silo {
  void optlist_free (DBoptlist *);

  struct OptVal {
    OptVal( int opt_id );
    std::variant<int, float, double> val;

    template < typename T >
    OptVal& operator=( T t ) noexcept {
      // val = t changes type in val to T, which is NOT what we want
      std::visit([=](auto& v) {v = t;}, val);
      return *this;
    }
  };

  struct OptList : public std::unordered_map<int,OptVal> {
  private:
    std::unique_ptr<DBoptlist, void(*)(DBoptlist*)> _p{nullptr, optlist_free};

  public:
    operator DBoptlist*();

    OptVal& operator[] ( int opt_id ) {
      emplace(opt_id, OptVal(opt_id));
      return at(opt_id);
    }
  };


}

#endif
