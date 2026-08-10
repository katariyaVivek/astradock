# Lesson 001: Vectors in Spacecraft Engineering

## Learning objective

A vector represents a quantity with both magnitude and direction. This lesson
explains how three-component vectors describe spacecraft physics, why their
coordinate frame and units are part of their meaning, how the principal vector
operations work, and how AstraDock verifies its first mathematical primitive.

## 1. Physical problem

A spacecraft has a position relative to an origin, moves with a velocity, and
accelerates when forces act on it. A rigid spacecraft also experiences torques.
Each of these quantities points somewhere, so one scalar is insufficient.

In three-dimensional Euclidean space a vector can be represented by three
components:

$$
\mathbf{v} = \begin{bmatrix}v_x & v_y & v_z\end{bmatrix}^{T}.
$$

The three numbers are meaningful only after the basis axes, coordinate frame,
and units are known.

## 2. Why vectors are needed

Spacecraft equations connect directional quantities:

- position $\mathbf{r}$ locates a vehicle relative to a chosen origin,
- velocity $\mathbf{v}$ describes the rate and direction of position change,
- acceleration $\mathbf{a}$ describes the rate of velocity change,
- force $\mathbf{F}$ changes translational motion,
- torque $\boldsymbol{\tau}$ changes rotational motion.

Using one `Vector3` arithmetic implementation avoids repeating component-wise
math in every later subsystem. It does not remove the engineering obligation to
label the frame and units.

## 3. Coordinate frames and frame dependence

Components are not the physical vector itself; they are the vector expressed in
a basis. The same physical velocity has different components in an
Earth-centred inertial (ECI) frame and a spacecraft body frame. Adding vectors
expressed in different frames is invalid even though each has three numbers.

Example variable names make the contract visible:

```cpp
Vector3 position_eci_m;
Vector3 velocity_eci_m_per_s;
Vector3 force_body_n;
```

`position_eci_m + force_body_n` is physically meaningless because both the
dimension and frame differ. The current `Vector3` type cannot detect this error,
so containing APIs and documentation must. Explicit frame types may be added
only when their benefit and conventions are understood.

## 4. Units

A vector's components share one physical dimension and unit. AstraDock uses SI
units internally:

| Quantity | Example symbol | SI unit |
| --- | --- | --- |
| Position | $\mathbf{r}$ | metre (m) |
| Velocity | $\mathbf{v}$ | metre per second (m/s) |
| Acceleration | $\mathbf{a}$ | metre per second squared (m/s^2) |
| Force | $\mathbf{F}$ | newton (N) |
| Torque | $\boldsymbol{\tau}$ | newton metre (N m) |
| Angular velocity | $\boldsymbol{\omega}$ | radian per second (rad/s) |

Vector arithmetic propagates units. Multiplying velocity by a time interval
produces a displacement; a dot product of two metre vectors has units of square
metres. Normalization divides out the magnitude and therefore produces a
dimensionless direction vector.

## 5. State variables, inputs, and outputs

At this foundation milestone the mathematical state is the component triplet
`(x, y, z)`, stored as double-precision values.

- Inputs: two vectors for binary operations, or one vector and a scalar for
  scaling.
- Outputs: a vector for addition, subtraction, scaling, cross product, and
  normalization; a scalar for dot product, squared norm, and norm.
- Invalid inputs: a zero divisor and the exactly zero vector when a direction is
  requested through normalization.

The vector type contains no time state, frame transform, measurement noise, or
spacecraft dynamics.

## 6. Magnitude and direction

The squared Euclidean norm is

$$
\lVert\mathbf{v}\rVert^2 = \mathbf{v}\cdot\mathbf{v}
= v_x^2 + v_y^2 + v_z^2.
$$

The magnitude is

$$
\lVert\mathbf{v}\rVert = \sqrt{v_x^2 + v_y^2 + v_z^2}.
$$

For example, $[3,4,0]$ has magnitude $5$. A nonzero vector's unit direction is

$$
\hat{\mathbf{v}} = \frac{\mathbf{v}}{\lVert\mathbf{v}\rVert}.
$$

The zero vector has no direction. AstraDock therefore throws `std::domain_error`
when asked to normalize it; returning another zero vector would hide an invalid
physical operation.

## 7. Addition, subtraction, and scalar operations

Addition and subtraction act component by component:

$$
\mathbf{a} + \mathbf{b}
= [a_x+b_x,\ a_y+b_y,\ a_z+b_z].
$$

These operations require compatible units and a common frame. Scalar
multiplication changes magnitude and may reverse direction when the scalar is
negative. Scalar division is multiplication by the reciprocal and is undefined
for a zero divisor.

## 8. Dot product

The dot product returns a scalar:

$$
\mathbf{a}\cdot\mathbf{b}
= a_xb_x + a_yb_y + a_zb_z
= \lVert\mathbf{a}\rVert\lVert\mathbf{b}\rVert\cos\theta.
$$

It measures alignment. A positive result means the vectors generally point in
the same direction, a negative result means they generally oppose one another,
and zero means nonzero vectors are perpendicular. Later, dot products will help
compute projections, directional errors, work, and geometric tests.

The input vectors must be expressed in the same orthonormal frame. The result's
unit is the product of the input units.

## 9. Cross product

The cross product returns a vector perpendicular to both inputs:

$$
\mathbf{a}\times\mathbf{b} =
\begin{bmatrix}
a_yb_z-a_zb_y \\
a_zb_x-a_xb_z \\
a_xb_y-a_yb_x
\end{bmatrix}.
$$

Its magnitude is

$$
\lVert\mathbf{a}\times\mathbf{b}\rVert
= \lVert\mathbf{a}\rVert\lVert\mathbf{b}\rVert\sin\theta.
$$

Direction follows the right-hand rule, so order matters:
$\mathbf{a}\times\mathbf{b}=-(\mathbf{b}\times\mathbf{a})$. In a right-handed
frame, $\hat{\mathbf{x}}\times\hat{\mathbf{y}}=\hat{\mathbf{z}}$.

Spacecraft examples include torque $\boldsymbol{\tau}=\mathbf{r}\times\mathbf{F}$
and angular momentum. Both inputs must use the same right-handed coordinate
frame, and the output uses the product of their units.

## 10. Governing assumptions and numerical method

The current implementation assumes:

- three-dimensional Euclidean geometry,
- a Cartesian, orthonormal, right-handed basis for geometric interpretations,
- double-precision IEEE 754 arithmetic,
- components already expressed in compatible frames and units.

Operations use their direct analytical component formulas. The norm uses the
standard library's three-argument `hypot`, which is generally safer against
intermediate overflow and underflow than directly evaluating the squared sum
before taking a square root.

Floating-point arithmetic is not exact. Tests use analytical exact cases where
appropriate and tolerance-based comparison for results affected by rounding.
The helper combines absolute tolerance near zero with relative tolerance as
values grow.

## 11. Known failure modes and common mistakes

- **Missing frame labels:** adding an ECI vector to a body-frame vector.
- **Mixed units:** adding metres to kilometres or velocity to position.
- **Wrong cross-product order:** reversing operands reverses the result.
- **Left-handed frame assumption:** the usual cross-product direction assumes a
  right-handed basis.
- **Normalizing zero:** zero magnitude cannot define a direction and throws.
- **Normalizing a tiny vector:** a nonzero but poorly scaled vector can amplify
  numerical noise; the owning physical algorithm must set an appropriate
  threshold.
- **Exact floating-point equality:** calculated values should normally be
  compared using a justified tolerance.
- **Overflow or non-finite inputs:** very large components, infinity, or NaN can
  produce non-finite results. Higher-level physical state validation will reject
  such inputs where its valid range is known.
- **Ignoring units of dot and cross products:** these operations do not
  automatically produce dimensionless quantities.

## 12. Verification strategy

The deterministic unit tests check:

1. known component-wise addition and subtraction,
2. scalar multiplication from either side and scalar division,
3. a dot product with a hand-calculated result,
4. right-hand-rule cross products of basis vectors,
5. squared norm and norm using a 3-4-5 triangle,
6. normalization components and unit magnitude,
7. the invariants $(\mathbf{a}\times\mathbf{b})\cdot\mathbf{a}=0$ and
   $(\mathbf{a}\times\mathbf{b})\cdot\mathbf{b}=0$,
8. explicit exceptions for zero division and zero-vector normalization.

These tests verify mathematical behavior without depending on orbital models or
random data.

## What you should now understand

1. Why do a vector's three components have no complete physical meaning without
   a coordinate frame and units?
2. Why is normalization undefined for the zero vector?
3. What geometric information does the sign of a dot product provide?
4. Why does reversing cross-product operands reverse the result?
5. Why must two vectors be in the same frame before adding them or taking their
   dot product?
6. Which spacecraft quantities are naturally represented as vectors?
7. Why are analytical cases and invariants both useful in unit tests?
