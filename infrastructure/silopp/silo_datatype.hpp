#ifndef _SILO_DATATYPE_
#define _SILO_DATATYPE_

namespace silo {
  template <typename Type>
  constexpr int datatype(Type* = nullptr) noexcept {
    using T = std::remove_const_t<Type>;
    if constexpr ( std::is_same_v<T, char> ) return DB_CHAR;
    else if ( std::is_same_v<T, short> ) return DB_SHORT;
    else if ( std::is_same_v<T, int> ) return DB_INT;
    else if ( std::is_same_v<T, long> ) return DB_LONG;
    else if ( std::is_same_v<T, long long> ) return DB_LONG_LONG;
    else if ( std::is_same_v<T, float> ) return DB_FLOAT;
    else if ( std::is_same_v<T, double> ) return DB_DOUBLE;
    else return DB_NOTYPE;
  }

  template <typename T>
  constexpr int datatype(const T&) noexcept {
    return datatype((T*)0);
  }
}

#endif
