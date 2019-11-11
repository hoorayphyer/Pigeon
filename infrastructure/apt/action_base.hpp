#ifndef  _APT_ACTION_BASE_HPP_
#define  _APT_ACTION_BASE_HPP_

#include "apt/index.hpp"
#include "apt/range.hpp"
#include <string>

namespace apt {
  template < int DGrid >
  struct ActionBase : public array<Range,DGrid> {
  private:
    std::string _name = "Unknown";

  public:
    ActionBase& setName(std::string name) { _name = name; return *this; }

    const auto& name() const noexcept { return _name; }

    virtual ~ActionBase() {};

    virtual ActionBase* Clone() const { return new ActionBase(*this); }; // covariant return types, see Modern C++ Design
  };
}

#endif
