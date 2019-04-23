#ifndef  _SILO_OPTLIST_HPP_
#define  _SILO_OPTLIST_HPP_

#include <variant>
#include <unordered_map>
#include <vector>
#include <memory>

typedef struct DBoptlist_ DBoptlist;
using RawList = DBoptlist;

namespace silo {
  void optlist_free (RawList *);

  struct OptList;

  struct OptVal {
  private:
    std::variant<int, float, double > _numeric;
    std::vector<int> _array; // NOTE variant is not allowed to hold reference, array, or void

  public:
    friend struct OptList;
    OptVal( int opt_id );

    // TODO does this change the type of val? Seems no from reference, which IS what we want
    template < typename T >
    std::enable_if_t< !std::is_same_v<T, std::vector<int>>, OptVal&> operator=( const T& t ) noexcept {
      std::visit([&t](auto& v) {v = t;}, _numeric);
      return *this;
    }

    OptVal& operator= ( const std::vector<int>& t ) {
      _array = t;
      return *this;
    }
  };

  struct OptList : public std::unordered_map<int,OptVal> {
  private:
    // used only to auto manage resources when converting to RawList*. No need to be copied in copy constructor
    std::unique_ptr<RawList, void(*)(RawList*)> _p{nullptr, optlist_free};

  public:
    operator RawList*();

    OptVal& operator[] ( int opt_id ) {
      emplace(opt_id, OptVal(opt_id));
      return at(opt_id);
    }

    OptList() = default;
    OptList( const OptList& other ) // simply call the copy constructor of unordered_map
      : std::unordered_map<int,OptVal>( static_cast<const std::unordered_map<int,OptVal>&>(other) ) {}
    OptList( OptList&& ) = default;
    ~OptList() = default;
  };


}

#endif
