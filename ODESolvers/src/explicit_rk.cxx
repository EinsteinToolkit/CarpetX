#include "explicit_rk.hxx"

#include <array>
#include <cstddef>
#include <numeric>
#include <limits>
#include <stdexcept>
#include <vector>

namespace ODESolvers {
namespace {

using I = std::int64_t;
constexpr RationalCoefficient Q(const I numerator,
                                const I denominator = 1) {
  return {numerator, denominator};
}

const ExplicitRKTableau rk4_tableau{
    ExplicitRKMethod::rk4,
    4,
    {{}, {Q(1, 2)}, {Q(0), Q(1, 2)}, {Q(0), Q(0), Q(1)}},
    {Q(1, 6), Q(1, 3), Q(1, 3), Q(1, 6)},
    {Q(0), Q(1, 2), Q(1, 2), Q(1)}};

const ExplicitRKTableau rkf78_tableau{
    ExplicitRKMethod::rkf78,
    7,
    {{},
     {Q(2, 27)},
     {Q(1, 36), Q(3, 36)},
     {Q(1, 24), Q(0), Q(3, 24)},
     {Q(20, 48), Q(0), Q(-75, 48), Q(75, 48)},
     {Q(1, 20), Q(0), Q(0), Q(5, 20), Q(4, 20)},
     {Q(-25, 108), Q(0), Q(0), Q(125, 108), Q(-260, 108), Q(250, 108)},
     {Q(31, 300), Q(0), Q(0), Q(0), Q(61, 225), Q(-2, 9), Q(13, 900)},
     {Q(2), Q(0), Q(0), Q(-53, 6), Q(704, 45), Q(-107, 9), Q(67, 90),
      Q(3)},
     {Q(-91, 108), Q(0), Q(0), Q(23, 108), Q(-976, 135), Q(311, 54),
      Q(-19, 60), Q(17, 6), Q(-1, 12)},
     {Q(2383, 4100), Q(0), Q(0), Q(-341, 164), Q(4496, 1025), Q(-301, 82),
      Q(2133, 4100), Q(45, 82), Q(45, 164), Q(18, 41)}},
    {Q(41, 840), Q(0), Q(0), Q(0), Q(0), Q(34, 105), Q(9, 35), Q(9, 35),
     Q(9, 280), Q(9, 280), Q(41, 840)},
    {Q(0), Q(2, 27), Q(1, 9), Q(1, 6), Q(5, 12), Q(1, 2), Q(5, 6), Q(1, 6),
     Q(2, 3), Q(1, 3), Q(1)}};

const ExplicitRKTableau dp87_tableau{
    ExplicitRKMethod::dp87,
    8,
    {{},
     {Q(1, 18)},
     {Q(1, 48), Q(1, 16)},
     {Q(1, 32), Q(0), Q(3, 32)},
     {Q(5, 16), Q(0), Q(-75, 64), Q(75, 64)},
     {Q(3, 80), Q(0), Q(0), Q(3, 16), Q(3, 20)},
     {Q(29443841, 614563906), Q(0), Q(0), Q(77736538, 692538347),
      Q(-28693883, 1125000000), Q(23124283, 1800000000)},
     {Q(16016141, 946692911), Q(0), Q(0), Q(61564180, 158732637),
      Q(22789713, 633445777), Q(545815736, 2771057229),
      Q(-180193667, 1043307555)},
     {Q(39632708, 573591083), Q(0), Q(0), Q(-433636366, 683701615),
      Q(-421739975, 2616292301), Q(100302831, 723423059),
      Q(790204164, 839813087), Q(800635310, 3783071287)},
     {Q(246121993, 1340847787), Q(0), Q(0), Q(-37695042795, 15268766246),
      Q(-309121744, 1061227803), Q(-12992083, 490766935),
      Q(6005943493, 2108947869), Q(393006217, 1396673457),
      Q(123872331, 1001029789)},
     {Q(-1028468189, 846180014), Q(0), Q(0), Q(8478235783, 508512852),
      Q(1311729495, 1432422823), Q(-10304129995, 1701304382),
      Q(-48777925059, 3047939560), Q(15336726248, 1032824649),
      Q(-45442868181, 3398467696), Q(3065993473, 597172653)},
     {Q(185892177, 718116043), Q(0), Q(0), Q(-3185094517, 667107341),
      Q(-477755414, 1098053517), Q(-703635378, 230739211),
      Q(5731566787, 1027545527), Q(5232866602, 850066563),
      Q(-4093664535, 808688257), Q(3962137247, 1805957418),
      Q(65686358, 487910083)},
     {Q(403863854, 491063109), Q(0), Q(0), Q(-5068492393, 434740067),
      Q(-411421997, 543043805), Q(652783627, 914296604),
      Q(11173962825, 925320556), Q(-13158990841, 6184727034),
      Q(3936647629, 1978049680), Q(-160528059, 685178525),
      Q(248638103, 1413531060), Q(0)}},
    {Q(14005451, 335480064), Q(0), Q(0), Q(0), Q(0),
     Q(-59238493, 1068277825), Q(181606767, 758867731),
     Q(561292985, 797845732), Q(-1041891430, 1371343529),
     Q(760417239, 1151165299), Q(118820643, 751138087),
     Q(-528747749, 2220607170), Q(1, 4)},
    {}};

void validate_coefficient(const RationalCoefficient coefficient) {
  if (coefficient.denominator <= 0)
    throw std::invalid_argument(
        "explicit RK coefficient denominator is not positive");
  if (coefficient.numerator == std::numeric_limits<I>::min())
    throw std::invalid_argument(
        "explicit RK coefficient numerator has unsupported minimum value");
}

bool same_rational(const RationalCoefficient left,
                   const RationalCoefficient right) {
  validate_coefficient(left);
  validate_coefficient(right);
  const auto left_divisor = std::gcd(left.numerator, left.denominator);
  const auto right_divisor = std::gcd(right.numerator, right.denominator);
  return left.numerator / left_divisor == right.numerator / right_divisor &&
         left.denominator / left_divisor ==
             right.denominator / right_divisor;
}

RationalCoefficient add_exact(const RationalCoefficient left,
                              const RationalCoefficient right) {
  validate_coefficient(left);
  validate_coefficient(right);
  const auto denominator_gcd =
      std::gcd(left.denominator, right.denominator);
  const auto left_multiplier = right.denominator / denominator_gcd;
  const auto right_multiplier = left.denominator / denominator_gcd;
  const auto checked_product = [](const I value, const I multiplier) {
    if (value != 0 &&
        (value > std::numeric_limits<I>::max() / multiplier ||
         value < std::numeric_limits<I>::min() / multiplier))
      throw std::overflow_error("exact RK rational validation overflow");
    return value * multiplier;
  };
  const I left_numerator = checked_product(left.numerator, left_multiplier);
  const I right_numerator = checked_product(right.numerator, right_multiplier);
  if ((right_numerator > 0 &&
       left_numerator > std::numeric_limits<I>::max() - right_numerator) ||
      (right_numerator < 0 &&
       left_numerator < std::numeric_limits<I>::min() - right_numerator))
    throw std::overflow_error("exact RK rational validation overflow");
  const I numerator = left_numerator + right_numerator;
  if (numerator == std::numeric_limits<I>::min())
    throw std::overflow_error("exact RK rational validation overflow");
  const I denominator = checked_product(left.denominator, left_multiplier);
  const auto divisor = std::gcd(numerator, denominator);
  return {numerator / divisor, denominator / divisor};
}

I checked_signed_product(const I left, const I right) {
  constexpr I minimum = std::numeric_limits<I>::min();
  constexpr I maximum = std::numeric_limits<I>::max();
  if (left > 0) {
    if ((right > 0 && left > maximum / right) ||
        (right < 0 && right < minimum / left))
      throw std::overflow_error("exact RK rational multiplication overflow");
  } else if (left < 0) {
    if ((right > 0 && left < minimum / right) ||
        (right < 0 && left < maximum / right))
      throw std::overflow_error("exact RK rational multiplication overflow");
  }
  return left * right;
}

RationalCoefficient multiply_exact(RationalCoefficient left,
                                   RationalCoefficient right) {
  validate_coefficient(left);
  validate_coefficient(right);
  const auto left_cancel = std::gcd(
      left.numerator < 0 ? -left.numerator : left.numerator,
      right.denominator);
  const auto right_cancel = std::gcd(
      right.numerator < 0 ? -right.numerator : right.numerator,
      left.denominator);
  left.numerator /= left_cancel;
  right.denominator /= left_cancel;
  right.numerator /= right_cancel;
  left.denominator /= right_cancel;
  return {checked_signed_product(left.numerator, right.numerator),
          checked_signed_product(left.denominator, right.denominator)};
}

RationalCoefficient sum_exact(
    const std::vector<RationalCoefficient> &coefficients) {
  RationalCoefficient result{0, 1};
  for (const auto coefficient : coefficients)
    result = add_exact(result, coefficient);
  return result;
}

void validate_against_canonical(const ExplicitRKTableau &tableau,
                                const ExplicitRKTableau &canonical) {
  for (std::size_t stage = 0; stage < tableau.a.size(); ++stage)
    for (std::size_t i = 0; i < tableau.a[stage].size(); ++i)
      if (tableau.a[stage][i].numerator !=
              canonical.a[stage][i].numerator ||
          tableau.a[stage][i].denominator !=
              canonical.a[stage][i].denominator)
        throw std::invalid_argument("explicit RK A coefficient mismatch");
  for (std::size_t i = 0; i < tableau.b.size(); ++i)
    if (tableau.b[i].numerator != canonical.b[i].numerator ||
        tableau.b[i].denominator != canonical.b[i].denominator)
      throw std::invalid_argument("explicit RK b coefficient mismatch");
  for (std::size_t i = 0; i < tableau.c.size(); ++i)
    if (tableau.c[i].numerator != canonical.c[i].numerator ||
        tableau.c[i].denominator != canonical.c[i].denominator)
      throw std::invalid_argument("explicit RK c coefficient mismatch");
}

void append_u8(std::vector<std::uint8_t> &bytes, const std::uint8_t value) {
  bytes.push_back(value);
}

void append_u32(std::vector<std::uint8_t> &bytes, const std::uint32_t value) {
  for (int shift = 24; shift >= 0; shift -= 8)
    bytes.push_back(static_cast<std::uint8_t>(value >> shift));
}

void append_u64(std::vector<std::uint8_t> &bytes, const std::uint64_t value) {
  for (int shift = 56; shift >= 0; shift -= 8)
    bytes.push_back(static_cast<std::uint8_t>(value >> shift));
}

void append_i64(std::vector<std::uint8_t> &bytes, const std::int64_t value) {
  append_u64(bytes, static_cast<std::uint64_t>(value));
}

void append_rational(std::vector<std::uint8_t> &bytes,
                     const RationalCoefficient coefficient) {
  validate_coefficient(coefficient);
  append_i64(bytes, coefficient.numerator);
  append_u64(bytes, static_cast<std::uint64_t>(coefficient.denominator));
}

std::vector<std::uint8_t>
serialize_tableau(const ExplicitRKTableau &tableau) {
  validate_explicit_rk_tableau(tableau);
  static constexpr std::array<std::uint8_t, 23> tag{{
      'C', 'X', '-', 'E', 'X', 'P', 'L', 'I', 'C', 'I', 'T', '-',
      'R', 'K', '-', 'T', 'A', 'B', 'L', 'E', 'A', 'U', '\0'}};
  std::vector<std::uint8_t> bytes;
  bytes.reserve(64U + tableau.a.size() * tableau.a.size() * 8U);
  bytes.insert(bytes.end(), tag.begin(), tag.end());
  append_u32(bytes, 1U);
  switch (tableau.method) {
  case ExplicitRKMethod::rk4:
    append_u8(bytes, 1U);
    break;
  case ExplicitRKMethod::rkf78:
    append_u8(bytes, 2U);
    break;
  case ExplicitRKMethod::dp87:
    append_u8(bytes, 3U);
    break;
  default:
    throw std::invalid_argument("unsupported explicit RK method");
  }
  append_u32(bytes, static_cast<std::uint32_t>(tableau.endpoint_order));
  append_u32(bytes, static_cast<std::uint32_t>(tableau.a.size()));
  append_u8(bytes, tableau.method == ExplicitRKMethod::dp87 ? 1U : 0U);
  append_u32(bytes, static_cast<std::uint32_t>(tableau.a.size()));
  for (const auto &row : tableau.a) {
    append_u32(bytes, static_cast<std::uint32_t>(row.size()));
    for (const auto coefficient : row)
      append_rational(bytes, coefficient);
  }
  append_u32(bytes, static_cast<std::uint32_t>(tableau.b.size()));
  for (const auto coefficient : tableau.b)
    append_rational(bytes, coefficient);
  append_u32(bytes, static_cast<std::uint32_t>(tableau.c.size()));
  for (const auto coefficient : tableau.c)
    append_rational(bytes, coefficient);
  return bytes;
}

constexpr std::uint32_t rotate_right(const std::uint32_t value,
                                     const int count) noexcept {
  return (value >> count) | (value << (32 - count));
}

ExplicitRKTableauFingerprint
sha256(const std::vector<std::uint8_t> &message) {
  static constexpr std::array<std::uint32_t, 64> round_constants{{
      0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U,
      0x3956c25bU, 0x59f111f1U, 0x923f82a4U, 0xab1c5ed5U,
      0xd807aa98U, 0x12835b01U, 0x243185beU, 0x550c7dc3U,
      0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U, 0xc19bf174U,
      0xe49b69c1U, 0xefbe4786U, 0x0fc19dc6U, 0x240ca1ccU,
      0x2de92c6fU, 0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU,
      0x983e5152U, 0xa831c66dU, 0xb00327c8U, 0xbf597fc7U,
      0xc6e00bf3U, 0xd5a79147U, 0x06ca6351U, 0x14292967U,
      0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU, 0x53380d13U,
      0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U,
      0xa2bfe8a1U, 0xa81a664bU, 0xc24b8b70U, 0xc76c51a3U,
      0xd192e819U, 0xd6990624U, 0xf40e3585U, 0x106aa070U,
      0x19a4c116U, 0x1e376c08U, 0x2748774cU, 0x34b0bcb5U,
      0x391c0cb3U, 0x4ed8aa4aU, 0x5b9cca4fU, 0x682e6ff3U,
      0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U,
      0x90befffaU, 0xa4506cebU, 0xbef9a3f7U, 0xc67178f2U}};
  std::array<std::uint32_t, 8> hash{{
      0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U, 0xa54ff53aU,
      0x510e527fU, 0x9b05688cU, 0x1f83d9abU, 0x5be0cd19U}};

  std::vector<std::uint8_t> padded(message);
  const std::uint64_t bit_count =
      static_cast<std::uint64_t>(message.size()) * 8U;
  padded.push_back(0x80U);
  while (padded.size() % 64U != 56U)
    padded.push_back(0U);
  append_u64(padded, bit_count);

  for (std::size_t block = 0; block < padded.size(); block += 64U) {
    std::array<std::uint32_t, 64> words{};
    for (std::size_t i = 0; i < 16U; ++i) {
      const std::size_t offset = block + 4U * i;
      words[i] = (static_cast<std::uint32_t>(padded[offset]) << 24) |
                 (static_cast<std::uint32_t>(padded[offset + 1U]) << 16) |
                 (static_cast<std::uint32_t>(padded[offset + 2U]) << 8) |
                 static_cast<std::uint32_t>(padded[offset + 3U]);
    }
    for (std::size_t i = 16U; i < words.size(); ++i) {
      const std::uint32_t s0 = rotate_right(words[i - 15U], 7) ^
                               rotate_right(words[i - 15U], 18) ^
                               (words[i - 15U] >> 3);
      const std::uint32_t s1 = rotate_right(words[i - 2U], 17) ^
                               rotate_right(words[i - 2U], 19) ^
                               (words[i - 2U] >> 10);
      words[i] = words[i - 16U] + s0 + words[i - 7U] + s1;
    }

    std::uint32_t a = hash[0];
    std::uint32_t b = hash[1];
    std::uint32_t c = hash[2];
    std::uint32_t d = hash[3];
    std::uint32_t e = hash[4];
    std::uint32_t f = hash[5];
    std::uint32_t g = hash[6];
    std::uint32_t h = hash[7];
    for (std::size_t i = 0; i < words.size(); ++i) {
      const std::uint32_t upper_sigma1 =
          rotate_right(e, 6) ^ rotate_right(e, 11) ^ rotate_right(e, 25);
      const std::uint32_t choice = (e & f) ^ ((~e) & g);
      const std::uint32_t temp1 =
          h + upper_sigma1 + choice + round_constants[i] + words[i];
      const std::uint32_t upper_sigma0 =
          rotate_right(a, 2) ^ rotate_right(a, 13) ^ rotate_right(a, 22);
      const std::uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
      const std::uint32_t temp2 = upper_sigma0 + majority;
      h = g;
      g = f;
      f = e;
      e = d + temp1;
      d = c;
      c = b;
      b = a;
      a = temp1 + temp2;
    }
    hash[0] += a;
    hash[1] += b;
    hash[2] += c;
    hash[3] += d;
    hash[4] += e;
    hash[5] += f;
    hash[6] += g;
    hash[7] += h;
  }

  ExplicitRKTableauFingerprint result{};
  for (std::size_t i = 0; i < hash.size(); ++i)
    for (std::size_t byte = 0; byte < 4U; ++byte)
      result[4U * i + byte] = static_cast<std::uint8_t>(
          hash[i] >> static_cast<int>(24U - 8U * byte));
  return result;
}

} // namespace

const ExplicitRKTableau &explicit_rk_tableau(const ExplicitRKMethod method) {
  switch (method) {
  case ExplicitRKMethod::rk4:
    return rk4_tableau;
  case ExplicitRKMethod::rkf78:
    return rkf78_tableau;
  case ExplicitRKMethod::dp87:
    return dp87_tableau;
  }
  throw std::invalid_argument("unsupported explicit RK method");
}

RationalCoefficient
explicit_rk_stage_fraction(const ExplicitRKMethod method,
                           const int stage_index) {
  static constexpr std::array<RationalCoefficient, 4> rk4_fractions{{
      Q(0), Q(1, 2), Q(1, 2), Q(1)}};
  static constexpr std::array<RationalCoefficient, 11> rkf78_fractions{{
      Q(0), Q(2, 27), Q(1, 9), Q(1, 6), Q(5, 12), Q(1, 2),
      Q(5, 6), Q(1, 6), Q(2, 3), Q(1, 3), Q(1)}};
  static constexpr std::array<RationalCoefficient, 13> dp87_fractions{{
      Q(0), Q(1, 18), Q(1, 12), Q(1, 8), Q(5, 16), Q(3, 8),
      Q(59, 400), Q(93, 200), Q(5490023248LL, 9719169821LL),
      Q(13, 20), Q(1201146811LL, 1299019798LL), Q(1), Q(1)}};
  if (stage_index <= 0)
    throw std::out_of_range("explicit RK stage index is out of range");
  const auto index = static_cast<std::size_t>(stage_index - 1);
  switch (method) {
  case ExplicitRKMethod::rk4:
    if (index < rk4_fractions.size())
      return rk4_fractions[index];
    break;
  case ExplicitRKMethod::rkf78:
    if (index < rkf78_fractions.size())
      return rkf78_fractions[index];
    break;
  case ExplicitRKMethod::dp87:
    if (index < dp87_fractions.size())
      return dp87_fractions[index];
    break;
  }
  throw std::out_of_range("explicit RK stage index is out of range");
}

ExplicitRKStagePoint
explicit_rk_stage_point(const ExplicitRKMethod method,
                        const ExplicitRKAdvanceFrame &frame,
                        const int stage_index) {
  if (frame.kind != ExplicitRKStageKind::primary &&
      frame.kind != ExplicitRKStageKind::fractional)
    throw std::invalid_argument(
        "explicit RK advance kind cannot be an endpoint probe");
  validate_coefficient(frame.begin_fraction);
  validate_coefficient(frame.extent_fraction);
  if (frame.begin_fraction.numerator < 0 ||
      frame.extent_fraction.numerator <= 0)
    throw std::invalid_argument("explicit RK advance fraction is invalid");
  const auto endpoint =
      add_exact(frame.begin_fraction, frame.extent_fraction);
  if (endpoint.numerator > endpoint.denominator)
    throw std::out_of_range(
        "explicit RK advance fraction exceeds the parent step");
  const auto &tableau = explicit_rk_tableau(method);
  const int stage_count = static_cast<int>(tableau.a.size());
  const auto local_fraction = explicit_rk_stage_fraction(method, stage_index);
  const auto parent_fraction =
      add_exact(frame.begin_fraction,
                multiply_exact(frame.extent_fraction, local_fraction));
  if (parent_fraction.numerator < 0 ||
      parent_fraction.numerator > parent_fraction.denominator)
    throw std::out_of_range(
        "explicit RK stage fraction exceeds the parent step");
  return {frame.kind, stage_index, stage_count, parent_fraction};
}

void validate_explicit_rk_tableau(const ExplicitRKTableau &tableau) {
  int expected_stages = 0;
  int expected_order = 0;
  const ExplicitRKTableau *canonical = nullptr;
  switch (tableau.method) {
  case ExplicitRKMethod::rk4:
    expected_stages = 4;
    expected_order = 4;
    canonical = &rk4_tableau;
    break;
  case ExplicitRKMethod::rkf78:
    expected_stages = 11;
    expected_order = 7;
    canonical = &rkf78_tableau;
    break;
  case ExplicitRKMethod::dp87:
    expected_stages = 13;
    expected_order = 8;
    canonical = &dp87_tableau;
    break;
  default:
    throw std::invalid_argument("unsupported explicit RK method");
  }
  if (tableau.endpoint_order != expected_order ||
      static_cast<int>(tableau.a.size()) != expected_stages ||
      static_cast<int>(tableau.b.size()) != expected_stages)
    throw std::invalid_argument("explicit RK descriptor shape/order mismatch");
  for (std::size_t stage = 0; stage < tableau.a.size(); ++stage) {
    if (tableau.a[stage].size() != stage)
      throw std::invalid_argument("explicit RK A row length mismatch");
    for (const auto coefficient : tableau.a[stage])
      validate_coefficient(coefficient);
  }
  for (const auto coefficient : tableau.b)
    validate_coefficient(coefficient);

  if (tableau.method == ExplicitRKMethod::dp87) {
    if (!tableau.c.empty())
      throw std::invalid_argument("DP87 must derive c from A rows");
  } else {
    if (tableau.c.size() != tableau.a.size())
      throw std::invalid_argument("explicit RK c descriptor length mismatch");
    for (const auto coefficient : tableau.c)
      validate_coefficient(coefficient);
  }
  if (tableau.method != ExplicitRKMethod::dp87 &&
      !same_rational(sum_exact(tableau.b), Q(1)))
    throw std::invalid_argument("exact explicit RK b sum is not one");
  if (tableau.method == ExplicitRKMethod::rkf78)
    for (std::size_t stage = 0; stage < tableau.a.size(); ++stage)
      if (!same_rational(sum_exact(tableau.a[stage]), tableau.c[stage]))
        throw std::invalid_argument("exact RKF78 A row sum does not equal c");
  // Execution reads the caller's exact integer pairs for the RK4 legacy
  // arithmetic. Require the audited canonical spelling, not merely a
  // mathematically equivalent fraction, before entering the kernel.
  validate_against_canonical(tableau, *canonical);
}

ExplicitRKTableauFingerprint
explicit_rk_tableau_fingerprint(const ExplicitRKMethod method) {
  return sha256(serialize_tableau(explicit_rk_tableau(method)));
}

} // namespace ODESolvers
