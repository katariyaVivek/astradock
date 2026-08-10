# Lesson 002: Central Two-Body Gravity

## Learning objective

This lesson develops the first physical model in AstraDock: the instantaneous
gravitational acceleration of a spacecraft near a spherical central body. It
connects Newton's law of gravitation to a frame-explicit vector equation,
explains the gravitational parameter `mu`, introduces the Cartesian orbital
state, and verifies the model with analytical Earth examples. It does not
propagate an orbit through time.

## 1. Physical problem

A spacecraft away from thrust and other disturbances still accelerates because
Earth attracts it. The spacecraft is not force-free in orbit; it is continually
falling toward Earth while its sideways velocity carries it around the planet.

For this milestone, the engineering question is deliberately narrow:

> Given the spacecraft's position relative to the central body's center, what
> is its instantaneous gravitational acceleration?

The answer becomes the dynamics input for a later numerical integrator, but no
integrator or propagation loop belongs in M02.

## 2. Position, velocity, acceleration, and orbital state

The position vector

$$
\mathbf r = [r_x, r_y, r_z]^T
$$

points from the central body's center of mass to the spacecraft. Its magnitude

$$
r = \lVert\mathbf r\rVert
$$

is the distance between their centers, not altitude above the surface.

The velocity vector

$$
\mathbf v = \frac{d\mathbf r}{dt}
$$

describes how position changes in the same frame. Position and velocity together
form a Cartesian orbital state:

$$
\mathbf x =
\begin{bmatrix}
\mathbf r \\
\mathbf v
\end{bmatrix}.
$$

Acceleration is the rate of change of velocity:

$$
\mathbf a = \frac{d\mathbf v}{dt} = \frac{d^2\mathbf r}{dt^2}.
$$

The M02 production function consumes only $\mathbf r$ and returns the
instantaneous $\mathbf a$. It does not store or advance $\mathbf x$.

## 3. Coordinate frame

Both position and acceleration are expressed in one
**central-body-centered inertial frame**:

- the origin is the central body's center of mass,
- the axes do not rotate with the central body,
- position and acceleration components use the same axes.

For an Earth scenario this is conceptually an Earth-centered inertial (ECI)
state. M02 does not define a particular epoch or realize an operational ECI
frame, and it does not transform between ECI and Earth-centered Earth-fixed
(ECEF) coordinates. Those transformations are a later milestone.

The current `Vector3` type cannot encode a frame, so the API parameter name
`position_central_body_inertial_m` and its documentation carry this contract.

## 4. Units

AstraDock uses SI units internally:

| Quantity | Symbol | Unit |
| --- | --- | --- |
| Position and radius | $\mathbf r$, $r$ | metre (m) |
| Velocity | $\mathbf v$ | metre per second (m/s) |
| Acceleration | $\mathbf a$ | metre per second squared (m/s^2) |
| Universal gravitational constant | $G$ | m^3/(kg s^2) |
| Central-body mass | $M$ | kilogram (kg) |
| Gravitational parameter | $\mu=GM$ | m^3/s^2 |

Supplying kilometres to an API that expects metres introduces errors by powers
of one thousand. Units therefore appear in constant and parameter names.

## 5. Newton's law of gravitation

For central-body mass $M$, spacecraft mass $m$, and center-to-center distance
$r$, Newton's law gives the force magnitude

$$
F = G\frac{Mm}{r^2}.
$$

The force points from the spacecraft toward the central body. Because
$\mathbf r/r$ points outward from the body to the spacecraft, the force vector
is

$$
\mathbf F = -G\frac{Mm}{r^2}\frac{\mathbf r}{r}.
$$

Dividing by spacecraft mass gives acceleration:

$$
\mathbf a = \frac{\mathbf F}{m}
= -GM\frac{\mathbf r}{r^3}.
$$

The spacecraft mass cancels. This is why, in the ideal model, objects at the
same position receive the same gravitational acceleration regardless of their
mass.

## 6. The gravitational parameter `mu`

Define

$$
\mu = GM.
$$

Then the model becomes

$$
\boxed{\mathbf a(\mathbf r)=-\mu\frac{\mathbf r}{\lVert\mathbf r\rVert^3}}.
$$

Aerospace software normally uses $\mu$ directly because spacecraft tracking
observes the dynamical product $GM$ much more accurately than separately
measuring the universal constant $G$ and the body's mass $M$. It also makes the
equations and units clearer.

For the introductory Earth examples, AstraDock defines the WGS 84 conventional
values:

```text
Earth mu                 = 3.986004418e14 m^3/s^2
Earth reference radius   = 6,378,137.0 m
```

The radius is the WGS 84 semi-major axis. M02 uses it as a documented spherical
reference radius for analytical examples; it is not a claim that Earth is
perfectly spherical or that every surface point is that distance from Earth's
center.

## Why is the denominator |r|^3 instead of |r|^2?

The inverse-square law specifies only the acceleration magnitude:

$$
\lVert\mathbf a\rVert = \frac{\mu}{r^2}.
$$

The vector equation also needs a direction. The outward radial unit vector is

$$
\hat{\mathbf r}=\frac{\mathbf r}{r}.
$$

Gravity points inward, so multiply the inverse-square magnitude by
$-\hat{\mathbf r}$:

$$
\mathbf a
= -\frac{\mu}{r^2}\frac{\mathbf r}{r}
= -\mu\frac{\mathbf r}{r^3}.
$$

One factor of $r$ normalizes the position vector; the other two factors create
the inverse-square magnitude. The vector equation still produces exactly
$\mu/r^2$ as the magnitude.

## 7. Why acceleration points toward the central body

The position vector points from the origin to the spacecraft. Multiplication by
a negative scalar reverses it, so

$$
\mathbf r\cdot\mathbf a < 0
$$

for every valid position. Position and acceleration are anti-parallel, giving

$$
\mathbf r\times\mathbf a=\mathbf 0.
$$

These two relationships are useful verification invariants for positions that
are not aligned with a coordinate axis.

## 8. Two-body assumptions

The name "two-body" refers to the central body and spacecraft. This initial
model assumes:

- the central body is spherical,
- gravity is represented by a point mass at the central body's center,
- spacecraft mass is negligible relative to central-body mass,
- the reference frame is inertial,
- there is no atmosphere,
- there is no aerodynamic drag,
- there is no third-body gravity from the Moon, Sun, or other objects,
- there is no oblateness or $J_2$ gravity term,
- there is no solar radiation pressure,
- there is no thrust or other non-gravitational acceleration.

No higher-fidelity correction is introduced in M02. The purpose is to understand
and verify the foundational central-gravity model first.

## 9. Inputs, outputs, and invalid states

The production API is conceptually:

```cpp
Vector3 two_body_acceleration(
    const Vector3& position_central_body_inertial_m,
    double gravitational_parameter_m3_per_s2);
```

- Input state used: central-body-relative inertial position in metres.
- Model parameter: positive, finite $\mu$ in m^3/s^2.
- Output: acceleration in the same inertial axes in m/s^2.

The origin $\mathbf r=\mathbf 0$ is a mathematical singularity: both direction
and inverse-square magnitude are undefined. AstraDock throws `std::domain_error`
rather than returning NaN or infinity. Non-finite position components and
non-positive or non-finite $\mu$ are also rejected. If a valid but extreme input
would produce a non-representable double-precision result, the function throws
`std::overflow_error`.

## 10. Numerical evaluation

This milestone evaluates the closed-form expression directly using
double-precision arithmetic. `Vector3::norm()` uses three-argument `hypot` to
reduce avoidable intermediate overflow and underflow. The implementation forms
$1/r$ once and multiplies it three times to obtain the scalar $-\mu/r^3$, then
scales the position vector.

There is no timestep, integration error, propagation convergence, or numerical
trajectory stability to analyze yet. Those topics begin in M03.

Known numerical and interface failure modes include:

- zero or extremely small radius near the singularity,
- non-finite inputs or results,
- negative or zero gravitational parameter,
- mixing metres with kilometres,
- using altitude where center-to-center radius is required,
- mixing coordinate frames,
- treating a central-body-fixed rotating vector as inertial,
- expecting surface-level or long-duration high-fidelity accuracy from a
  point-mass model.

## 11. Circular orbital speed

For a circular orbit, velocity is tangent to the orbit and perpendicular to the
radius. The required inward centripetal acceleration magnitude is

$$
a_c=\frac{v^2}{r}.
$$

Equating it to two-body gravity,

$$
\frac{v^2}{r}=\frac{\mu}{r^2},
$$

gives

$$
\boxed{v_\text{circular}=\sqrt{\frac{\mu}{r}}}.
$$

This is an analytical condition, not a propagation algorithm. M02 keeps the
calculation in its lesson and tests rather than adding another production API.

## 12. Specific orbital energy concept

The ideal two-body specific mechanical energy is energy per unit spacecraft
mass:

$$
\varepsilon=\frac{\lVert\mathbf v\rVert^2}{2}-\frac{\mu}{r}.
$$

Its units are m^2/s^2, equivalent to joules per kilogram. Negative energy
describes a bound orbit, zero is the ideal escape boundary, and positive energy
describes an unbound trajectory. In an ideal propagated two-body orbit,
$\varepsilon$ should remain constant. AstraDock will implement and monitor this
in later propagation/invariant milestones, not in M02.

## 13. Specific angular momentum concept

Specific angular momentum is

$$
\mathbf h=\mathbf r\times\mathbf v,
$$

with units m^2/s. Its direction is normal to the orbital plane according to the
right-hand rule. Central gravity produces no torque about the origin because
$\mathbf r\times\mathbf a=\mathbf 0$, so $\mathbf h$ is conserved in the ideal
two-body problem. Production helpers and conservation checks belong with orbit
propagation and invariants in M04/M05.

## 14. Worked 500 km Low Earth Orbit example

Use the WGS 84 reference radius and an altitude of 500 km:

```text
Earth reference radius  R_E = 6,378,137 m
altitude                 h   =   500,000 m
orbital radius           r   = 6,878,137 m
Earth mu                 mu  = 3.986004418e14 m^3/s^2
```

The two-body gravitational acceleration magnitude is

$$
\frac{\mu}{r^2}
= \frac{3.986004418\times10^{14}}{(6.878137\times10^6)^2}
\approx 8.426\ \text{m/s}^2.
$$

The analytical circular speed is

$$
\sqrt{\frac{\mu}{r}}
= \sqrt{\frac{3.986004418\times10^{14}}{6.878137\times10^6}}
\approx 7,612.6\ \text{m/s}
\approx 7.613\ \text{km/s}.
$$

Gravity at 500 km is still strong. Orbital flight is sustained free fall, not an
absence of gravity.

## 15. Earth-surface sanity check

At the reference radius, the point-mass result is

$$
\frac{\mu}{R_E^2}\approx9.7983\ \text{m/s}^2,
$$

which is close to the familiar 9.8 m/s^2. This is a sanity check, not a
high-fidelity surface-gravity prediction. Real apparent gravity varies with
latitude, elevation, Earth's rotation, ellipsoidal shape, and nonuniform mass
distribution.

## 16. Verification strategy

Deterministic analytical tests verify:

1. a position on positive X produces acceleration on negative X,
2. acceleration magnitude equals $\mu/r^2$,
3. X, Y, and Z axis cases are equivalent and point toward the origin,
4. an arbitrary position satisfies $\mathbf r\times\mathbf a\approx0$ and
   $\mathbf r\cdot\mathbf a<0$,
5. doubling radius reduces acceleration magnitude by a factor of four,
6. Earth reference-radius acceleration is approximately 9.8 m/s^2,
7. the 500 km analytical circular speed is approximately 7.613 km/s,
8. zero position and invalid physical parameters fail explicitly.

These checks compare directly with the governing physics. They do not depend on
a trajectory looking visually plausible.

## 17. Model limitations

M02 computes one instantaneous acceleration vector. It does not provide:

- Euler or RK4 integration,
- an orbital propagation loop,
- ECI/ECEF or LVLH transforms,
- Earth rotation,
- $J_2$ or higher-order gravity,
- atmospheric drag,
- third-body forces,
- attitude or rigid-body dynamics,
- guidance, navigation, control, sensors, or estimation.

## Sources for the Earth constants

- The [NGA WGS 84 reference page](https://earth-info.nga.mil/index.php?action=wgs84&dir=wgs84)
  identifies the 6,378,137.0 m semi-major axis and links the defining standard.
- NASA's [6-DOF simulation check-case appendix](https://nescacademy.nasa.gov/src/flightsim/Reports/NASA-TM-2015-218675-EOM_checkcase_appendices.pdf)
  tabulates the WGS 84 gravitational parameter and reference radii in SI units.

## What you should now understand

1. Why does the acceleration vector point opposite the position vector?
2. Why does the vector equation contain $r^3$ even though gravity follows an
   inverse-square magnitude law?
3. What does $\mu=GM$ represent, and why is it used directly?
4. Why is altitude not the same quantity as orbital radius?
5. Why are both position and velocity required to define an orbital state?
6. Why is circular velocity perpendicular to the radius vector?
7. What do specific energy and specific angular momentum reveal about an ideal
   orbit?
8. Which physical effects are absent from the two-body model?
