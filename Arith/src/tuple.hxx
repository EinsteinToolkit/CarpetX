#ifndef CARPETX_ARITH_TUPLE_HXX
#define CARPETX_ARITH_TUPLE_HXX

#include <cstddef>
#include <type_traits>
#include <utility>

#ifndef __CUDACC__

#include <tuple>

template <typename... Args> using std_tuple = std::tuple<Args...>;

template <typename... Args>
constexpr decltype(auto) std_make_tuple(Args &&...args) {
  return std::make_tuple(std::forward<Args>(args)...);
}

template <std::size_t I, typename... Args>
using std_tuple_element = std::tuple_element<I, Args...>;

template <typename... Args> using std_tuple_size = std::tuple_size<Args...>;

template <typename... Args>
inline constexpr bool std_tuple_size_v = std::tuple_size_v<Args...>;

#else

#include <thrust/tuple.h>

template <typename... Args> using std_tuple = thrust::tuple<Args...>;

template <typename... Args>
__device__ __host__ constexpr decltype(auto) std_make_tuple(Args &&...args) {
  return thrust::make_tuple(std::forward<Args>(args)...);
}

template <std::size_t I, typename... Args>
using std_tuple_element = thrust::tuple_element<I, Args...>;

template <typename... Args> using std_tuple_size = thrust::tuple_size<Args...>;

template <typename... Args>
inline constexpr bool std_tuple_size_v = std_tuple_size<Args...>::value;

namespace std {
template <std::size_t I, typename... Args>
constexpr decltype(auto) get(const std_tuple<Args...> &x) {
  return thrust::get<I>(x);
}
template <std::size_t I, typename... Args>
constexpr decltype(auto) get(std_tuple<Args...> &&x) {
  return thrust::get<I>(std::move(x));
}
} // namespace std

#endif

#endif // #ifndef CARPETX_ARITH_TUPLE_HXX
