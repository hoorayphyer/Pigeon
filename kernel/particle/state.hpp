#pragma once

#include <bitset>
#include <cstdint>
#include <tuple>
#include <type_traits>

#include "apt/bit_manip.hpp"
#include "apt/ternary.hpp"
#include "particle/flags_predef.hpp"
#include "particle/species_predef.hpp"

namespace particle {
template <int N, typename T = std::uint32_t>
struct proxy_int {
  using proxy_int_t = T;
  static constexpr int nbits_v = N;
  static_assert(N <= sizeof(T) * 8);

  constexpr proxy_int() = default;
  constexpr proxy_int(proxy_int_t a) noexcept : _data(a) {}
  constexpr operator proxy_int_t() const noexcept { return _data; }

 protected:
  proxy_int_t _data{};
};

struct flagbits : public std::bitset<16> {
 private:
  using base_type = std::bitset<16>;

 public:
  using proxy_int_t = std::uint32_t;
  static constexpr int nbits_v = 16;

  using base_type::bitset;

  // explicit to avoid ambiguity in operator==(int)
  explicit operator proxy_int_t() const noexcept {
    return static_cast<proxy_int_t>(to_ulong());
  }

  constexpr const auto& operator[](flag a) const {
    return static_cast<const base_type*>(this)->operator[](
        static_cast<std::underlying_type_t<flag>>(a));
  }

  decltype(auto) operator[](flag a) {
    return static_cast<base_type*>(this)->operator[](
        static_cast<std::underlying_type_t<flag>>(a));
  }
};

struct migrcode : proxy_int<5> {
  using proxy_int::proxy_int;

  template <int D>
  constexpr migrcode& decode() noexcept {
    // fill all leading bits with the value from sign bit. (two's compliment)
    apt::setbits<nbits_v - 1, sizeof(_data) * 8 - nbits_v + 1>(
        _data, ~(_data >> (nbits_v - 1)) + 1);
    _data += (apt::pow3(D) - 1) / 2;
    return *this;
  }

  template <int D>
  constexpr void encode() noexcept {
    // NOTE due to two's compliment, no need to explicitly set the sign bit at
    // nbits_v - 1 location
    _data = static_cast<proxy_int_t>(_data - (apt::pow3(D) - 1) / 2);
    // still explicitly reset all leading bits, just to be consistent with the
    // result read off from a state object
    _data &= ((1u << nbits_v) - 1);
  }
};

using freebits = proxy_int<0>;

using pid = proxy_int<1 + 12 + 23, std::uint64_t>;  // particle id

template <typename T>
struct to_int {
  using type = typename T::proxy_int_t;
};

template <>
struct to_int<species> {
  using type = std::underlying_type_t<species>;
};

template <typename T>
using to_int_t = typename to_int<T>::type;

template <typename T>
struct nbits {
  static constexpr int value = T::nbits_v;
};

template <>
struct nbits<species> {
  static constexpr int value = 3;
};

template <typename T>
constexpr int nbits_v = nbits<T>::value;
}  // namespace particle

namespace particle {
struct layout {
 private:
  using ordering = std::tuple<species, migrcode, flagbits, freebits, pid>;
  static constexpr auto sizing =
      std::make_tuple(nbits_v<species>, nbits_v<migrcode>, nbits_v<flagbits>,
                      nbits_v<freebits>, nbits_v<pid>);

 public:
  template <typename Attr, std::size_t I = 0>
  static constexpr std::size_t size() noexcept {
    static_assert((I < std::tuple_size_v<ordering>), "unknown Attr");
    using elm_t = std::tuple_element_t<I, ordering>;
    if constexpr (std::is_same_v<Attr, elm_t>)
      return std::get<I>(sizing);
    else
      return size<Attr, I + 1>();
  }

  template <typename Attr, std::size_t I = 0>
  static constexpr std::size_t begin() noexcept {
    static_assert((I < std::tuple_size_v<ordering>), "unknown Attr");
    using elm_t = std::tuple_element_t<I, ordering>;
    if constexpr (std::is_same_v<Attr, elm_t>)
      return 0;
    else
      return size<elm_t>() + begin<Attr, I + 1>();
  }
};

template <typename Ptc, typename T>
struct StateExpression {
 private:
  template <typename Attr>
  inline void set_impl(const Attr& attr) noexcept {
    apt::setbits<layout::begin<Attr>(), layout::size<Attr>()>(
        state(), static_cast<to_int_t<Attr>>(attr));
  }

  template <typename Attr>
  inline Attr get_impl() const noexcept {
    return static_cast<Attr>(
        apt::getbits<layout::begin<Attr>(), layout::size<Attr>()>(state()));
  }

 public:
  using self_t = StateExpression<Ptc, T>;

  constexpr T& state() noexcept { return static_cast<Ptc&>(*this).state(); }
  constexpr T state() const noexcept {
    return static_cast<const Ptc&>(*this).state();
  }

  template <typename Attr, typename... Attrs>
  inline std::enable_if_t<!std::is_same_v<Attr, migrcode>, self_t&> set(
      const Attr& attr, const Attrs&... attrs) noexcept {
    static_assert(std::is_same_v<Attr, species> or
                  std::is_same_v<Attr, flagbits> or
                  std::is_same_v<Attr, freebits> or std::is_same_v<Attr, pid> or
                  std::is_same_v<Attr, flag>);

    if constexpr (!std::is_same_v<Attr, flag>)
      set_impl(attr);
    else {
      auto x = static_cast<std::underlying_type_t<flag>>(attr);
      state() |=
          (static_cast<std::size_t>(1) << (x + layout::begin<flagbits>()));
    }

    if constexpr (sizeof...(Attrs) > 0) set(attrs...);
    return *this;
  }

  template <typename Migr_t, int D>
  inline std::enable_if_t<std::is_same_v<Migr_t, migrcode>, self_t&> set(
      int x) noexcept {
    migrcode m(x);
    m.encode<D>();
    set_impl(m);
    return *this;
  }

  inline self_t& reset(const flag& fl) noexcept {
    auto x = static_cast<std::underlying_type_t<flag>>(fl);
    // NOTE must cast 1 to std::size_t to avoid being wiped out
    state() &=
        ~(static_cast<std::size_t>(1) << (x + layout::begin<flagbits>()));
    return *this;
  }

  template <typename Attr>
  inline self_t& reset() noexcept {
    set_impl(Attr{});
    return *this;
  }

  template <typename Attr>
  inline std::enable_if_t<!std::is_same_v<Attr, migrcode>, Attr> get()
      const noexcept {
    return get_impl<Attr>();
  }

  template <typename Migr_t, int D>
  inline std::enable_if_t<std::is_same_v<Migr_t, migrcode>, to_int_t<migrcode>>
  get() const noexcept {
    migrcode m(get_impl<migrcode>());
    return m.decode<D>();
  }

  inline bool is(const flag& fl) const noexcept {
    auto x = static_cast<std::underlying_type_t<flag>>(fl);
    return (state() >> layout::begin<flagbits>()) &
           (static_cast<std::size_t>(1) << x);
  }

  template <flag F>
  inline self_t& assign(bool b) noexcept {
    return b ? set(F) : reset(F);
  }
};

}  // namespace particle
