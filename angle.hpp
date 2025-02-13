#pragma once
#include "pi.hpp"
#include <cmath>

namespace math {
template <typename T, usize Dimensions> struct Direction;
namespace detail {
template <typename T, usize Dimensions> struct VecBase;
}

enum AngleOrder { XYZ, /* XZY, YXZ, YZX, ZXY, ZYX unsupported for now */ };

template <typename T> struct SinCosPair {
  T sine, cosine;
};

// Implements convenience functions for both of these
#define IMPL_ANGLE                                                             \
  operator T &() { return value; }                                             \
  operator T const &() const { return value; }                                 \
  T sine() const { return std::sin(torad(value)); }                            \
  T cosine() const { return std::cos(torad(value)); }                          \
  SinCosPair<T> sincos() const { return {sine(), cosine()}; }

template <typename T> struct Radians {
  using Number = T;
  static inline constexpr T torad(T const &in) { return in; }
  static inline constexpr T todeg(T const &in) { return rad2deg(in); }
  Radians<T> squash() const { return std::fmod(value, (T)(2.0 * pi)); }
  T value;
  IMPL_ANGLE;
};

template <typename T> struct Degree {
  using Number = T;
  static inline constexpr T torad(T const &in) { return deg2rad(in); }
  static inline constexpr T todeg(T const &in) { return in; }
  Radians<T> squash() const { return std::fmod(value, (T)360.0); }
  T value;
  IMPL_ANGLE;
};

template <typename AngleNumber, AngleOrder Order = XYZ> struct EulerAng {
  using Number = AngleNumber::Number;
  union {
    struct {
      AngleNumber pitch{}, yaw{}, roll{};
    };
    struct {
      AngleNumber x, y, z;
    };
    AngleNumber data[3];
  };

  // Conversion between Formats:
  EulerAng<Degree<Number>, Order> as_degrees() const {
    return {AngleNumber::todeg(pitch), AngleNumber::todeg(yaw),
            AngleNumber::todeg(roll)};
  }

  EulerAng<Radians<Number>, Order> as_radians() const {
    return {AngleNumber::torad(pitch), AngleNumber::torad(yaw),
            AngleNumber::torad(roll)};
  }

  // Conversion between mathematical concepts
  Direction<Number, 3> forward() const {
    auto const [sp, cp] = pitch.sincos();
    auto const [sy, cy] = yaw.sincos();

    return {cp * cy, cp * sy, -sp};
  }

  Direction<Number, 3>::Relatives relatives() const {
    auto const [sp, cp] = pitch.sincos();
    auto const [sy, cy] = yaw.sincos();
    auto const [sr, cr] = roll.sincos();

    return {
        {cp * cy, cp * sy, -sp},
        {-(sr * sp * cy) + (cr * sy), -(sr * sp * sy) + -(cr * cy), -(sr * cp)},
        {cr * sp * cy + -sr * -sy, cr * sp * sy + -sr * cy, cr * cp}};
  }

  detail::VecBase<Number, 4> to_quat() const {
    auto const [sp, cp] = pitch.sincos();
    auto const [sy, cy] = yaw.sincos();
    auto const [sr, cr] = roll.sincos();

    auto const sr_x_cp = sr * cp, cr_x_sp = cr * sp;
    auto const cr_x_cp = cr * cp, sr_x_sp = sr * sp;

    return {sr_x_cp * cy - cr_x_sp * sy, cr_x_sp * cy + sr_x_cp * sy,
            cr_x_cp * sy - sr_x_sp * cy, cr_x_cp * cy + sr_x_sp * sy};
  }
};
} // namespace math
