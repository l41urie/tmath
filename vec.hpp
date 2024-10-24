#pragma once
#include "meta.hpp"
#include <array>
#include <cmath>

namespace math {
/// Represents a Point in N-Dimensional space
template <typename T, usize Dimensions> struct Point;
/// Represents a Direction in N-Dimensional space.
template <typename T, usize Dimensions> struct Direction;
/// Represents the difference between 2 Points in N-Dimensional space
template <typename T, usize Dimensions> struct Delta;
/// Represents a 3d-plane
template <typename T> struct Plane3D;
template <typename T = f64>
inline Plane3D<T> plane(Point<T, 3> origin, const Direction<T, 3> &fwd,
                        const Direction<T, 3> &up);

namespace detail {
template <typename T, usize M, usize N> struct MatBase {
  using NumberType = T;
  static auto constexpr SizeM = M;
  static auto constexpr SizeN = N;

  std::array<std::array<T, SizeN>, SizeM> data;

  operator T const *() const { return (T const *)data.data(); }
  operator T *() { return (T *)data.data(); }

#define BASE_OP(name, op, arg)                                                 \
  MatBase<T, M, N> name(arg) const {                                           \
    MatBase<T, M, N> result;                                                   \
    for (usize i = 0; i < M; ++i)                                              \
      for (usize j = 0; j < N; ++j)                                            \
        result.data[i][j] = this->data[i][j] op;                               \
    return result;                                                             \
  }

  BASE_OP(hadamard, *rhs.data, MatBase const &rhs);
  BASE_OP(sub, -rhs.data[i][j], MatBase const &rhs);
  BASE_OP(add, +rhs.data[i][j], MatBase const &rhs);

  BASE_OP(scale, *scalar, T const &scalar);
  BASE_OP(divide, / scalar, T const &scalar);

  T l1norm() const {
    T result{};
    for (usize i = 0; i < M; ++i)
      for (usize j = 0; j < N; ++j)
        result += std::abs(this->data[i][j]);
    return result;
  }

  // Given a 1xN matrix, this is equivalent to computing the dot product
  T frobenius_inner_product(MatBase<T, M, N> const &rhs) const {
    T result{};
    for (usize i = 0; i < M; ++i)
      for (usize j = 0; j < N; ++j)
        result += this->data[i][j] * rhs.data[i][j];
    return result;
  }

  // also "squared magnitude/length"
  T frobenius_norm_sq() const { return frobenius_inner_product(*this); }
};

template <typename T, usize Dimensions> struct VecData {
  MatBase<T, 1, Dimensions> data{};
};

template <typename T> struct VecData<T, 2> {
  union {
    MatBase<T, 1, 2> data{};
    struct {
      T x, y;
    };
  };
};

template <typename T> struct VecData<T, 3> {
  union {
    MatBase<T, 1, 3> data{};
    struct {
      T x, y, z;
    };
  };

protected:
#define IMPL_CROSS(in)                                                         \
  Vtd::Del cross(in const &rhs) const { return {VecBase<T, 3>::cross({rhs})}; }

  VecData<T, 3> cross(VecData<T, 3> const &rhs) const {
    return {this->y * rhs.z - this->z * rhs.y,
            this->x * rhs.z - this->z * rhs.x,
            this->x * rhs.y - this->y * rhs.x};
  }
};

template <typename T> struct VecData<T, 4> {
  union {
    MatBase<T, 1, 4> data{};
    struct {
      T x, y, z, w;
    };
  };
};

template <typename T, usize Dimensions>
struct VecBase : public VecData<T, Dimensions> {
  T &operator[](usize idx) { return this->data.data[0][idx]; }
  T const &operator[](usize idx) const { return this->data.data[0][idx]; }

protected:
  T dot(VecBase const &rhs) const {
    return this->data.frobenius_inner_product(rhs.data);
  }

  T len_sqr() const { return this->data.frobenius_norm_sq(); }
  T len() const { return std::sqrt(len_sqr()); }
  VecBase<T, Dimensions> norm() const { return {this->data.divide(len())}; }
};

// Deliberately left empty, anything that can be applied to _ALL_ dimensions
// directly goes into the struct
template <typename T, usize Dimensions>
struct PointBase : public VecBase<T, Dimensions> {};
template <typename T, usize Dimensions>
struct DirectionBase : public VecBase<T, Dimensions> {};
template <typename T, usize Dimensions>
struct DeltaBase : public VecBase<T, Dimensions> {};

template <typename T, usize Dimensions> struct VecTypeDesc {
  using P = Point<T, Dimensions>;
  using Dir = Direction<T, Dimensions>;
  using Del = Delta<T, Dimensions>;
};

// Specializations for some Dimensions that are commonly used.
template <typename T> struct PointBase<T, 3> : public VecBase<T, 3> {
  /// Projects this point onto a plane,
  /// located at |org|, oriented by |u| and |v|
  Point<T, 2> project(Point<T, 3> const &org, Direction<T, 3> const &u,
                      Direction<T, 3> const &v) const {
    auto d = (Point<T, 3>{*this} - org);
    return d.project(u, v);
  }

  /// Projects this point onto a plane,
  Point<T, 2> project(Plane3D<T> const &plane) const {
    return project(plane.location, plane.u, plane.v);
  }
};

template <typename T> struct DirectionBase<T, 3> : public VecBase<T, 3> {
  using Vtd = VecTypeDesc<T, 3>;
  IMPL_CROSS(Vtd::Dir);

  struct Relatives {
    Direction<T, 3> fwd, right, up;

    Plane3D<T> as_plane(Point<T, 3> const &location) const {
      return plane(location, fwd, up);
    }
  };

  Relatives relatives(bool pitch_up = true) const {
    if (fabs(this->x) < 1e-6 && fabs(this->y) < 1e-6) {
      // pitch 90 degrees up/down from identity
      if (pitch_up) {
        return {*this, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
      }
      return {*this, {0.0, -1.0, 0.0}, {-this->z, 0.0, 0.0}};
    }
    auto const right = Delta<T, 3>{this->y, this->x, 0.0}.norm();
    return {*this, right, this->cross(right).norm()};
  }
};

template <typename T> struct DeltaBase<T, 3> : public VecBase<T, 3> {
  using Vtd = VecTypeDesc<T, 3>;
  IMPL_CROSS(Vtd::Del);

  /// Projects a Delta onto a 2d plane oriented in U and V
  Point<T, 2> project(Direction<T, 3> const &u,
                      Direction<T, 3> const &v) const {
    return {this->dot(u), this->dot(v)};
  }
};
} // namespace detail

// clang-format off
#define OPERATOR(in, out, op, fn) out operator op (in const &rhs) const { return {this->data. fn (rhs.data)}; }
#define OPERATOR_ASSIGN(in, op, fn) void operator op (in const &rhs) { this->data = this->data. fn (rhs.data); }
#define OPERATOR_PRIM(in, out, op, fn) out operator op (in const &rhs) const { return {this->data. fn (rhs)}; }
#define OPERATOR_ASSIGN_PRIM(in, op, fn) void operator op (in const &rhs) { this->data = this->data. fn (rhs); }
#define IMPL_DOT(input) T dot(input const &rhs) const { return this->data.dot(rhs.data); }
// clang-format on

template <typename T, usize Dimensions>
struct Point : public detail::PointBase<T, Dimensions> {
  using Vtd = detail::VecTypeDesc<T, Dimensions>;

  /// Lowers the dimension by discarding the last scalar.
  Point<T, Dimensions - 1> flatten() {
    Point<T, Dimensions - 1> result;
    for (usize i = 0; i < Dimensions - 1; ++i)
      result[i] = (*this)[i];
    return result;
  }

  OPERATOR(Vtd::P, Vtd::Del, -, sub);
  OPERATOR(Vtd::Del, Vtd::P, -, sub);
  OPERATOR_ASSIGN(Vtd::Del, -=, sub);
  OPERATOR(Vtd::Del, Vtd::P, +, add);
  OPERATOR_ASSIGN(Vtd::Del, +=, add);

  /// Returns the Direction to another point.
  Vtd::Dir to(Vtd::P const &rhs) { return {(rhs - *this).norm()}; }
};

template <typename T, usize Dimensions>
struct Direction : public detail::DirectionBase<T, Dimensions> {
  using Vtd = detail::VecTypeDesc<T, Dimensions>;

  /// Lowers the dimension by discarding the last scalar.
  /// This correctly preserves The length
  Direction<T, Dimensions - 1> flatten() {
    Delta<T, Dimensions - 1> result;
    for (usize i = 0; i < Dimensions - 1; ++i)
      result[i] = (*this)[i];
    return result.norm();
  }

  OPERATOR_PRIM(T, Vtd::Del, *, scale);
  OPERATOR_PRIM(T, Vtd::Del, /, divide);

  IMPL_DOT(Vtd::Del);
  IMPL_DOT(Vtd::Dir);
};

template <typename T, usize Dimensions>
struct Delta : public detail::DeltaBase<T, Dimensions> {
  using Vtd = detail::VecTypeDesc<T, Dimensions>;

  /// Lowers the dimension by discarding the last scalar.
  /// The length changes here.
  Delta<T, Dimensions - 1> flatten() {
    Delta<T, Dimensions - 1> result;
    for (usize i = 0; i < Dimensions - 1; ++i)
      result[i] = (*this)[i];
    return result;
  }

  OPERATOR(Vtd::Del, Vtd::Del, +, add);
  OPERATOR(Vtd::P, Vtd::P, +, add);
  OPERATOR_ASSIGN(Vtd::Del, +=, add);

  OPERATOR_PRIM(T, Vtd::Del, *, scale);
  OPERATOR_ASSIGN_PRIM(T, *=, scale);
  OPERATOR_PRIM(T, Vtd::Del, /, divide);
  OPERATOR_ASSIGN_PRIM(T, /=, divide);

  using detail::VecBase<T, Dimensions>::len;
  using detail::VecBase<T, Dimensions>::len_sqr;

  OPERATOR(T, bool, <, shorter);
  OPERATOR(Vtd::Del, bool, <, shorter);
  OPERATOR(T, bool, >, longer);
  OPERATOR(Vtd::Del, bool, >, longer);

  IMPL_DOT(Vtd::Del);

  /// Normalizes the Delta to a Direction, length is guaranteed to be 1
  Vtd::Dir norm() const { return {detail::VecBase<T, Dimensions>::norm()}; }

  /// Forces the length to a |scalar|.
  Vtd::Del with_len(T const &scalar) const { return norm() * scalar; }

  /// Ensures the length isn't longer than |scalar|
  Vtd::Del clamped_to_len(T const &scalar) const {
    return *this < scalar ? *this : with_len(scalar);
  }
};

template <typename T> struct Plane3D {
  Direction<T, 3> u, v;
  Point<T, 3> location;
};

/// represents a point in Space
template <typename T = f64> using P2 = Point<T, 2>;
template <typename T = f64> using P3 = Point<T, 3>;
template <typename T = f64> using P4 = Point<T, 4>;

template <typename T = f64> using Del2 = Delta<T, 2>;
template <typename T = f64> using Del3 = Delta<T, 3>;
template <typename T = f64> using Del4 = Delta<T, 4>;

template <typename T = f64> math::P2<T> inline point(T const &x, T const &y) {
  return {x, y};
}

template <typename T = f64>
math::P3<T> inline point(T const &x, T const &y, T const &z) {
  return {x, y, z};
}

template <typename T = f64>
math::P4<T> inline point(T const &x, T const &y, T const &z, T const &w) {
  return {x, y, z, w};
}

/// Anything that represents something relative from a Point
template <typename T = f64> Del3<T> inline dir(T const &x, T const &y) {
  return {x, y};
}

template <typename T = f64>
Del3<T> inline dir(T const &x, T const &y, T const &z) {
  return {x, y, z};
}

template <typename T>
inline Plane3D<T> plane(Point<T, 3> location, Direction<T, 3> const &fwd,
                        Direction<T, 3> const &up) {
  auto u = fwd.cross(up).norm();

  // n = forward
  // const auto u = forward.Cross( Vector( 0, 0, 1 ) ).Normalize();
  // const auto v = u.Cross( n ).Normalize();
  return Plane3D<T>{
      .location = location,
      .u = u,
      .v = u.cross(up).norm(),
  };
}
} // namespace math
