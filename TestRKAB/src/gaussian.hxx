#ifndef TEST_RKAB_GAUSSIAN_HXX
#define TEST_RKAB_GAUSSIAN_HXX

#include <cctk.h>

#include <cmath>

namespace TestRKAB::gauss {

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE phi(T W, T A, T t, T x, T y,
                                             T z) noexcept -> T {
  using std::sqrt, std::cosh, std::sinh;

  const auto r{sqrt(x * x + y * y + z * z)};

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    return (2 * A * t) / (exp(pow(t, 2) / (2. * pow(W, 2))) * pow(W, 2));
  } else {
    return (A *
            (exp(-0.5 * pow(t - sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)), 2) /
                 pow(W, 2)) -
             exp(-0.5 * pow(t + sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)), 2) /
                 pow(W, 2)))) /
           sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  }
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Pi(T W, T A, T t, T x, T y,
                                            T z) noexcept -> T {
  using std::sqrt, std::cosh, std::sinh;

  const auto r{sqrt(x * x + y * y + z * z)};

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    return (2 * A * (-t + W) * (t + W)) /
           (exp(pow(t, 2) / (2. * pow(W, 2))) * pow(W, 4));
  } else {
    return (2 * A *
            (sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) *
                 cosh((t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                      pow(W, 2)) -
             t * sinh((t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                      pow(W, 2)))) /
           (exp((pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) /
                (2. * pow(W, 2))) *
            pow(W, 2) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
  }
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Dx(T W, T A, T t, T x, T y,
                                            T z) noexcept -> T {
  using std::sqrt, std::cosh, std::sinh;

  const auto r{sqrt(x * x + y * y + z * z)};

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    return 0;
  } else {
    return (A * x *
            (pow(W, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2) +
             t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) +
             exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                 pow(W, 2)) *
                 t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) -
             exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                 pow(W, 2)) *
                 (pow(W, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)))) /
           (exp((pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2) +
                 2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                (2. * pow(W, 2))) *
            pow(W, 2) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
  }
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Dy(T W, T A, T t, T x, T y,
                                            T z) noexcept -> T {
  using std::sqrt, std::cosh, std::sinh;

  const auto r{sqrt(x * x + y * y + z * z)};

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    return 0;
  } else {
    return (A * y *
            (pow(W, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2) +
             t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) +
             exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                 pow(W, 2)) *
                 t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) -
             exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                 pow(W, 2)) *
                 (pow(W, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)))) /
           (exp((pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2) +
                 2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                (2. * pow(W, 2))) *
            pow(W, 2) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
  }
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Dz(T W, T A, T t, T x, T y,
                                            T z) noexcept -> T {
  using std::sqrt, std::cosh, std::sinh;

  const auto r{sqrt(x * x + y * y + z * z)};

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    return 0;
  } else {
    return (A * z *
            (pow(W, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2) +
             t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) +
             exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                 pow(W, 2)) *
                 t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) -
             exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                 pow(W, 2)) *
                 (pow(W, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)))) /
           (exp((pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2) +
                 2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                (2. * pow(W, 2))) *
            pow(W, 2) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
  }
}

} // namespace TestRKAB::gauss

#endif // TEST_RKAB_GAUSSIAN_HXX
