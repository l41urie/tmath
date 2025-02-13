# tmath - A C++ Mathematics Library for Enhanced Type Safety

tmath is a C++20 maths library designed to provide easy-to-use abstractions for mathematical computations while enforcing stricter type safety than most of it's counterparts.
By leveraging C++'s type system, tmath minimizes runtime errors and enhances readability by using well-defined types for mathematical entities.

This library is particularly beneficial for applications in graphics, physics simulations, and any domain where precision and correctness in mathematical operations are paramount.

## Key Features

- **Type-Safe Mathematical Constructs**: Each mathematical entity (e.g., points, directions, deltas, angles) is represented by a distinct type, preventing unintended operations between incompatible types.
- **Intuitive API**: The library provides a High-level API that allows users without a strong mathematical background to write mathematical code confidently reducing the risk of errors.

## Example Usage

Here’s a quick example demonstrating the usage of tmath with annotated types:

```cpp
using namespace math;
auto const z = point(0.0, 0.0, 0.0); // Point<3>, a point in 3d-space
auto const up = point(0.0, 0.0, 100.0);

// Calculate the difference between two points
auto const diff = z - up; // Delta<3>, 3D Delta between 2 points

// Normalize the delta to get a direction
auto const norm = diff.norm(); // Dir<3>, a 3D Direction, internally, the length is guaranteed to be = 1 here.

// Convert the direction to Euler angles in degrees
auto const angle = norm.as_euler().as_degrees(); // EulerAng<Degrees>
```

### Preventing errors

To illustrate the library's error prevention capabilities, consider the following examples (in the context of the code above):
```cpp
z + diff; // Ok
z + up; // Errors
z + norm; // Errors
norm.norm(); // Errors

```

1. `z + diff` This operation is valid as it represents the addition of a delta vector to a point in three-dimensional space. In the context of affine geometry, this operation is well-defined and corresponds to the translation of the point z by the vector diff.

2. `z + up` This results in an error because the addition of two points is not a mathematically valid operation. In Euclidean space, points represent locations and do not possess an inherent additive structure. Thus, the operation lacks a meaningful interpretation within the framework of vector spaces.

3. `z + norm` This operation also yields an error, as it attempts to add a direction vector (represented by norm, which has a unit length) to a point. Conceptually, a direction is a vector that indicates orientation but does not possess a position in space. Therefore, the addition of a direction to a point is not a valid operation, as it conflates distinct mathematical entities.

4. `norm.norm()` This operation generates an error because normalizing a direction vector is redundant. By definition, a normalized vector has a magnitude of one, and thus applying the normalization operation again is nonsensical.
