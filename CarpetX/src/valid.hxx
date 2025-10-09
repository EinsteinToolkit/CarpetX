#ifndef CARPETX_CARPETX_VALID_HXX
#define CARPETX_CARPETX_VALID_HXX

#include <tuple.hxx>

#include <AMReX_Box.H>

#include <yaml-cpp/yaml.h>
#include <zlib.h>

#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

namespace CarpetX {

struct valid_t {
  bool valid_int, valid_outer, valid_ghosts;
  constexpr valid_t() : valid_t(false) {}

  explicit constexpr valid_t(bool b)
      : valid_int(b), valid_outer(b), valid_ghosts(b) {}
  friend constexpr valid_t operator~(const valid_t &x) {
    valid_t r;
    r.valid_int = !x.valid_int;
    r.valid_outer = !x.valid_outer;
    r.valid_ghosts = !x.valid_ghosts;
    return r;
  }
  friend constexpr valid_t operator&(const valid_t &x, const valid_t &y) {
    valid_t r;
    r.valid_int = x.valid_int && y.valid_int;
    r.valid_outer = x.valid_outer && y.valid_outer;
    r.valid_ghosts = x.valid_ghosts && y.valid_ghosts;
    return r;
  }
  friend constexpr valid_t operator|(const valid_t &x, const valid_t &y) {
    valid_t r;
    r.valid_int = x.valid_int || y.valid_int;
    r.valid_outer = x.valid_outer || y.valid_outer;
    r.valid_ghosts = x.valid_ghosts || y.valid_ghosts;
    return r;
  }
  valid_t &operator&=(const valid_t &x) { return *this = *this & x; }
  valid_t &operator|=(const valid_t &x) { return *this = *this | x; }

  constexpr bool valid_all() const {
    return valid_int && valid_outer && valid_ghosts;
  }
  constexpr bool valid_any() const {
    return valid_int || valid_outer || valid_ghosts;
  }

  friend bool operator==(const valid_t &x, const valid_t &y) {
    return std::make_tuple(x.valid_int, x.valid_outer, x.valid_ghosts) ==
           std::make_tuple(y.valid_int, y.valid_outer, y.valid_ghosts);
  }
  friend bool operator<(const valid_t &x, const valid_t &y) {
    return std::make_tuple(x.valid_int, x.valid_outer, x.valid_ghosts) <
           std::make_tuple(y.valid_int, y.valid_outer, y.valid_ghosts);
  }

  std::string explanation() const;

  friend std::ostream &operator<<(std::ostream &os, const valid_t v);
  operator std::string() const;
  friend YAML::Emitter &operator<<(YAML::Emitter &yaml, const valid_t v);
};

inline constexpr valid_t make_valid_int() {
  valid_t valid;
  valid.valid_int = true;
  return valid;
}
inline constexpr valid_t make_valid_outer() {
  valid_t valid;
  valid.valid_outer = true;
  return valid;
}
inline constexpr valid_t make_valid_ghosts() {
  valid_t valid;
  valid.valid_ghosts = true;
  return valid;
}
inline constexpr valid_t make_valid_all() { return ~valid_t(); }

} // namespace CarpetX
namespace std {
using namespace CarpetX;
template <> struct equal_to<valid_t> {
  constexpr bool operator()(const valid_t &x, const valid_t &y) const {
    return std::make_tuple(x.valid_int, x.valid_outer, x.valid_ghosts) ==
           std::make_tuple(y.valid_int, y.valid_outer, y.valid_ghosts);
  }
};
template <> struct less<valid_t> {
  constexpr bool operator()(const valid_t &x, const valid_t &y) const {
    return std::make_tuple(x.valid_int, x.valid_outer, x.valid_ghosts) <
           std::make_tuple(y.valid_int, y.valid_outer, y.valid_ghosts);
  }
};
} // namespace std
namespace CarpetX {

class why_valid_t {
  valid_t valid;
  std::function<std::string()> why_int, why_outer, why_ghosts;

public:
  // The constructor that doesn't give a reason should never be called
  why_valid_t() = delete;
  why_valid_t(const std::function<std::string()> &why)
      : why_valid_t(false, why) {}
  why_valid_t(bool b, const std::function<std::string()> &why)
      : why_valid_t(valid_t(b), why) {}
  why_valid_t(const valid_t &val, const std::function<std::string()> &why)
      : valid(val), why_int(why), why_outer(why), why_ghosts(why) {}

  const valid_t &get() const { return valid; }

  void set(const valid_t &which, const valid_t &val,
           const std::function<std::string()> &why) {
    valid = (valid & ~which) | (val & which);
    if (which.valid_int)
      why_int = why;
    if (which.valid_outer)
      why_outer = why;
    if (which.valid_ghosts)
      why_ghosts = why;
  }
  void set_all(const valid_t &val, const std::function<std::string()> &why) {
    set(valid_t(true), val, why);
  }
  void set_int(bool b, const std::function<std::string()> &why) {
    set(make_valid_int(), valid_t(b), why);
  }
  void set_outer(bool b, const std::function<std::string()> &why) {
    set(make_valid_outer(), valid_t(b), why);
  }
  void set_ghosts(bool b, const std::function<std::string()> &why) {
    set(make_valid_ghosts(), valid_t(b), why);
  }
  void set_invalid(const valid_t &which,
                   const std::function<std::string()> &why) {
    set(which, valid_t(false), why);
  }
  void set_valid(const valid_t &which,
                 const std::function<std::string()> &why) {
    set(which, valid_t(true), why);
  }

  std::string explanation() const;

  friend std::ostream &operator<<(std::ostream &os, const why_valid_t &why);
  operator std::string() const;
  friend YAML::Emitter &operator<<(YAML::Emitter &yaml, const why_valid_t &why);
};

////////////////////////////////////////////////////////////////////////////////

struct checksum_t {
  valid_t where;
  uLong crc;
  checksum_t() = default;
  inline checksum_t(const valid_t &where)
      : where(where), crc(crc32(0, nullptr, 0)) {}
  template <typename T> inline void add(const T &x) {
    crc = crc32(crc, static_cast<const Bytef *>(static_cast<const void *>(&x)),
                sizeof x);
  }

  friend bool operator==(const checksum_t &x, const checksum_t &y) {
    return x.where == y.where && x.crc == y.crc;
  }
  friend bool operator!=(const checksum_t &x, const checksum_t &y) {
    return !(x == y);
  }

  friend std::ostream &operator<<(std::ostream &os, const checksum_t &x) {
    return os << "checksum_t{where:" << x.where << ",crc:0x" << std::hex
              << std::setfill('0') << std::setw(8) << x.crc << "}";
  }
  operator std::string() const {
    std::ostringstream buf;
    buf << *this;
    return buf.str();
  }
};

struct tiletag_t {
  int patch, level, component;
  int gi, vi, tl;
  tiletag_t() = delete;

  friend bool operator==(const tiletag_t &x, const tiletag_t &y) {
    return std::make_tuple(x.patch, x.level, x.component, x.gi, x.vi, x.tl) ==
           std::make_tuple(y.patch, y.level, y.component, y.gi, y.vi, y.tl);
  }
  friend bool operator<(const tiletag_t &x, const tiletag_t &y) {
    return std::make_tuple(x.patch, x.level, x.component, x.gi, x.vi, x.tl) <
           std::make_tuple(y.patch, y.level, y.component, y.gi, y.vi, y.tl);
  }

  friend std::ostream &operator<<(std::ostream &os, const tiletag_t &x) {
    return os << "tiletag_t{"
              << "patch:" << x.patch << ","
              << "level:" << x.level << ","
              << "component:" << x.component << ","
              << "gi:" << x.gi << ","
              << "vi:" << x.vi << ","
              << "tl:" << x.tl << "}";
  }
  operator std::string() const {
    std::ostringstream buf;
    buf << *this;
    return buf.str();
  }
};

typedef std::map<tiletag_t, checksum_t> checksums_t;

checksums_t calculate_checksums(
    const std::vector<std::vector<std::vector<valid_t> > > &will_write);
void check_checksums(const checksums_t &checksums,
                     const std::function<std::string()> &where);

////////////////////////////////////////////////////////////////////////////////

template <typename T> struct ipoison_t;

template <> struct ipoison_t<CCTK_COMPLEX> {
#if defined CCTK_REAL_PRECISION_4
  const std::uint32_t value[2] = {0xffc00000UL + 0xdead, 0xffc00000UL + 0xdead};
#elif defined CCTK_REAL_PRECISION_8
  const std::uint64_t value[2] = {0xfff8000000000000ULL + 0xdeadbeef,
                                  0xfff8000000000000ULL + 0xdeadbeef};
#endif
  static_assert(sizeof value == sizeof(CCTK_COMPLEX));
};

template <> struct ipoison_t<CCTK_REAL> {
#if defined CCTK_REAL_PRECISION_4
  const std::uint32_t value[1] = {0xffc00000UL + 0xdead};
#elif defined CCTK_REAL_PRECISION_8
  const std::uint64_t value[1] = {0xfff8000000000000ULL + 0xdeadbeef};
#endif
  static_assert(sizeof value == sizeof(CCTK_REAL));
};

template <> struct ipoison_t<CCTK_INT> {
  const std::uint32_t value[1] = {0xdeadbeef};
  static_assert(sizeof value == sizeof(CCTK_INT));
};

template <typename T> class poison_value_t {
public:
  poison_value_t() = default;

  CCTK_HOST CCTK_DEVICE bool is_poison(const T &val) const {
    // no memcmp on CUDA :-(
    // return std::memcmp(&val, &ipoison.value, sizeof(val)) == 0;

    typedef decltype(ipoison.value) array_type;
    typedef typename std::remove_extent<array_type>::type const_type;
    const size_t array_size = std::extent<array_type>::value;

    const_type *ival = reinterpret_cast<const_type *>(&val);
    bool retval = false;
    for (size_t i = 0; i < array_size; ++i)
      if (ival[i] == ipoison.value[i])
        retval = true;
    return retval;
  }

  void set_to_poison(T &val) const {
    std::memcpy(static_cast<void *>(&val), &ipoison.value, sizeof(val));
  }

  void set_to_poison(void *val, size_t count) const {
    set_to_poison(reinterpret_cast<T *>(val), count);
  }

  void set_to_poison(T *val, size_t count) const {
    for (size_t i = 0; i < count; ++i)
      set_to_poison(val[i]);
  }

private:
  const ipoison_t<T> ipoison;
};

} // namespace CarpetX

#endif // #ifndef CARPETX_CARPETX_VALID_HXX
