#pragma once
#include "meta.hpp"

namespace math {
static constexpr f64 pi = 3.14159265358979323846;
static constexpr f32 pi32 = (f32)pi;

template <typename T> constexpr T deg2rad(T deg) {
  return deg * ((T)pi / (T)180.0);
}

template <typename T> constexpr T rad2deg(T rad) {
  return rad * ((T)180.0 / (T)pi);
}
} // namespace math