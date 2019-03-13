#ifndef _OLD_FIELDUPDATER_ADAPTER_HPP_
#define _OLD_FIELDUPDATER_ADAPTER_HPP_

#include "../parameters.hpp"

namespace field {
  template < typename, int, int > struct Field;
}

namespace ofsa {
  // only applicable in 2D log spherical
  struct OldFieldUpdater {
  private:
    void* _pfu = nullptr; // to be casted

  public:
    using field_type = field::Field<double,3,2>;
    OldFieldUpdater();
    ~OldFieldUpdater();

    // TODO factor of 4\pi on J ?
    void operator( field_type& E,
                   field_type& B,
                   const field_type& J,
                   field_type::element_t dt );
  };
}

#endif
