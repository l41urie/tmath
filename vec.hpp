#pragma once
#include "meta.hpp"
#include <cmath>
#include <numeric>

namespace math {
// PX Represents a point in 3D-space
// DirX Represents a Normalized Vector, length is guaranteed to be = 1
// DelX Represents a Delta between 2 vectors

template <typename T = f64> struct Vec4;
template <typename T = f64> struct P4;
template <typename T = f64> struct Dir4;
template <typename T = f64> struct Del4;

template <typename T = f64> struct Vec3;
template <typename T = f64> struct P3;
template <typename T = f64> struct Dir3;
template <typename T = f64> struct Del3;

template <typename T = f64> struct Vec2;
template <typename T = f64> struct P2;
template <typename T = f64> struct Dir2;
template <typename T = f64> struct Del2;

// Base struct.
template <typename T, usize SIZE> struct VecBase {
  using NumberType = T;
  static auto constexpr Size = SIZE;
  T data[Size]{};

  operator T const *() const { return data; }
  operator T *() { return data; }

#define BASE_OP(name, op, arg)                                                 \
  VecBase<T, SIZE> name(arg) const {                                           \
    VecBase<T, SIZE> result;                                                   \
    for (usize i = 0; i < SIZE; ++i)                                           \
      result[i] = this->data[i] op;                                            \
    return result;                                                             \
  }

  BASE_OP(hadamard, *rhs.data[i], VecBase const &rhs);
  BASE_OP(sub, -rhs.data[i], VecBase const &rhs);
  BASE_OP(add, +rhs.data[i], VecBase const &rhs);

  BASE_OP(scale, *scalar, T const &scalar);
  BASE_OP(divide, / scalar, T const &scalar);

  T l1norm() const { return std::accumulate(data, data + Size, (T)0); }

  T dot(VecBase const &rhs) const { return this->hadamard(rhs).l1norm(); }

  T len_sqr() const { return dot(*this); }
  T len() const { return std::sqrt(len_sqr()); }

  bool shorter(T const &scalar) const { return len_sqr() < (scalar * scalar); }
  bool shorter(VecBase const &rhs) const { return len_sqr() < rhs.len_sqr(); }

  VecBase<T, SIZE> norm() const { return divide(len()); }
};

// clang-format off
#define OPERATOR(in, out, op, fn) out operator op (in const &rhs) const { return {this->data. fn (rhs.data)}; }
#define OPERATOR_ASSIGN(in, op, fn) void operator op (in const &rhs) { this->data = this->data. fn (rhs.data); }

#define OPERATOR_PRIM(in, out, op, fn) out operator op (in const &rhs) const { return {this->data. fn (rhs)}; }
#define OPERATOR_ASSIGN_PRIM(in, op, fn) void operator op (in const &rhs) { this->data = this->data. fn (rhs); }

#define IMPLEMENT_LENGTH(del_type) \
  T len_sqr() const { return this->data.len_sqr(); } \
  T len() const { return this->data.len(); } \
  OPERATOR(T, bool, <, shorter) \
  OPERATOR(del_type, bool, <, shorter)

#define IMPL_DOT(input) T dot(input const &rhs) const { return this->data.dot(rhs.data); }

#define IMPL_NORM(out_dir, out_del) \
  out_dir norm() const { return {this->data.norm()}; } \
  out_del with_len(T const &scalar) const { return norm() * scalar; } \
  \
  Del3<T> clamped_to_len(T const &scalar) { \
    if (this->data.shorter(scalar)) \
      return *this; \
    return with_len(scalar); \
  }

#define OPERATORS_POINT(p_type, dir_type, del_type) \
  OPERATOR(p_type, del_type, -, sub); \
  OPERATOR(del_type, p_type, -, sub); \
  OPERATOR_ASSIGN(del_type, -=, sub); \
  OPERATOR(del_type, p_type, +, add); \
  OPERATOR_ASSIGN(del_type, +=, add); \
  dir_type to(p_type const &rhs) { return (rhs - *this).norm(); }

/*technically, division could be considered valid
but doing this is silly and most likely an error
*/
#define OPERATORS_DIR(p_type, dir_type, del_type) \
  OPERATOR_PRIM(T, del_type, *, scale); \
  IMPL_DOT(del_type); \
  IMPL_DOT(dir_type);


#define OPERATORS_DEL(p_type, dir_type, del_type)\
  OPERATOR(del_type, del_type, +, add); \
  OPERATOR(p_type, p_type, +, add); \
  OPERATOR_ASSIGN(del_type, +=, add); \
  OPERATOR_PRIM(T, del_type, *, scale); \
  OPERATOR_ASSIGN_PRIM(T, *=, scale); \
  OPERATOR_PRIM(T, del_type, /, divide); \
  OPERATOR_ASSIGN_PRIM(T, /=, divide); \
  IMPLEMENT_LENGTH(del_type); \
  IMPL_DOT(del_type); \
  IMPL_NORM(dir_type, del_type);

// clang-format on

// Definitions for 4d-space
template <typename T> struct Vec4 {
  using NumberType = T;
  union {
    VecBase<T, 4> data;
    T x, y, z, w;
  };
};

template <typename T> struct P4 : public Vec4<T> {
  OPERATORS_POINT(P4<T>, Dir4<T>, Del4<T>);
};

template <typename T> struct Dir4 : public Vec4<T> {
  OPERATORS_DIR(P4<T>, Dir4<T>, Del4<T>);
};

template <typename T> struct Del4 : public Vec4<T> {
  OPERATORS_DEL(P4<T>, Dir4<T>, Del4<T>);
};

// Definitions for 3d-space
template <typename T> struct Vec3 {
  using NumberType = T;
  union {
    VecBase<T, 3> data;
    T x, y, z;
  };

protected:
  Vec3<T> cross(Vec3 const &rhs) const {
    return {this->y * rhs.z - this->z * rhs.y,
            this->x * rhs.z - this->z * rhs.x,
            this->x * rhs.y - this->y * rhs.x};
  }
};

template <typename T> struct P3 : public Vec3<T> {
  OPERATORS_POINT(P3<T>, Dir3<T>, Del3<T>);

  P2<T> as2d() { return {this->x, this->y}; }
  // TODO: Projection stuff?
};

template <typename T> struct Dir3 : public Vec3<T> {
  OPERATORS_DIR(P3<T>, Dir3<T>, Del3<T>);

  Del3<T> cross(Del3<T> const &rhs) const { return {Vec3<T>::cross(rhs)}; }
  Del3<T> cross(Dir3<T> const &rhs) const { return {Vec3<T>::cross(rhs)}; }

  struct Relatives {
    Dir3<T> fwd, right, up;
  };

  Relatives relatives(bool pitch_up = true) const {
    if (fabs(this->x) < 1e-6 && fabs(this->y) < 1e-6) {
      // pitch 90 degrees up/down from identity
      if (pitch_up) {
        return {*this, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
      }
      return {*this, {0.0, -1.0, 0.0}, {-this->z, 0.0, 0.0}};
    }
    auto const right = Del3<T>{this->y, this->x, 0.0}.norm();
    return {*this, right, this->cross(right).norm()};
  }

  Dir2<T> as2d() { return Del2<T>{this->x, this->y}.norm(); }
};

template <typename T> struct Del3 : public Vec3<T> {
  OPERATORS_DEL(P3<T>, Dir3<T>, Del3<T>);

  Del3<T> cross(Del3<T> const &rhs) const { return Vec3<T>::cross(rhs); }
  Del2<T> as2d() { return {this->x, this->y}; }
};

// Definitions for 2d-space

template <typename T> struct Vec2 {
  using NumberType = T;
  union {
    VecBase<T, 2> data;
    T x, y;
  };
};

template <typename T> struct P2 : public Vec2<T> {
  OPERATORS_POINT(P2<T>, Dir2<T>, Del2<T>);
};

template <typename T> struct Dir2 : public Vec2<T> {
  OPERATORS_DIR(P2<T>, Dir2<T>, Del2<T>);
};

template <typename T> struct Del2 : public Vec2<T> {
  OPERATORS_DEL(P2<T>, Dir2<T>, Del2<T>);
};

// represents a point in Space
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

// Anything that represents something relative from a Point
template <typename T = f64> Del3<T> inline dir(T const &x, T const &y) {
  return {x, y};
}

template <typename T = f64>
Del3<T> inline dir(T const &x, T const &y, T const &z) {
  return {x, y, z};
}
} // namespace math