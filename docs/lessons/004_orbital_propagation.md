# Lesson 004: Two-Body Orbital Propagation

## Learning objective

M02 answered which gravitational acceleration acts at one position. M03
answered how a computer advances a state from its derivative. M04 connects
those ideas so a Cartesian spacecraft state moves through time, then uses
analytical circular-orbit geometry and conservation laws to judge whether the
numerical trajectory is trustworthy.

```text
gravity provides acceleration
integrator advances the state
repeated integration creates a trajectory
conservation laws diagnose numerical fidelity
```

## 1. Why orbital motion is second-order

Newton's second law determines acceleration. For the ideal two-body model,

\[
\ddot{\mathbf r}
=
-\mu\frac{\mathbf r}{\lVert\mathbf r\rVert^3}.
\]

The second derivative appears because gravity tells us how velocity changes,
while position is obtained only after velocity is also advanced. Position
alone is insufficient to predict the future: two spacecraft at the same point
but with different velocities follow different trajectories.

## 2. Converting to a first-order system

Euler and RK4 operate on first-order systems,

\[
\dot{\mathbf x}=f(t,\mathbf x).
\]

Define a six-component translational state:

\[
\mathbf x=
\begin{bmatrix}
\mathbf r\\
\mathbf v
\end{bmatrix}
=
\begin{bmatrix}
r_x&r_y&r_z&v_x&v_y&v_z
\end{bmatrix}^{T}.
\]

Its derivative is

\[
\dot{\mathbf x}
=
\begin{bmatrix}
\dot{\mathbf r}\\
\dot{\mathbf v}
\end{bmatrix}
=
\begin{bmatrix}
\mathbf v\\
-\mu\mathbf r/\lVert\mathbf r\rVert^3
\end{bmatrix}.
\]

The first relationship, \(\dot{\mathbf r}=\mathbf v\), says velocity is the
instantaneous rate of position change. The second,
\(\dot{\mathbf v}=\mathbf a\), says acceleration is the instantaneous rate of
velocity change.

AstraDock's derivative returns the same algebraic `CartesianState` shape used
for the physical state. In a physical state its slots contain position in
metres and velocity in m/s. In a derivative, the first slot contains position
rate in m/s and the second contains velocity rate in m/s^2. This explicit
contract avoids a large dimensional-type framework while keeping the RK4
arithmetic readable.

## 3. Frame and units

M04 uses an idealized **Earth-centered inertial Cartesian frame**:

- the origin is Earth's center of mass,
- axes are treated as non-rotating for the ideal model,
- position, velocity, acceleration, and angular momentum components use those
  same axes.

This is a simplified inertial frame, not a realized operational ECI convention
with an epoch. There is no Earth-fixed frame, Earth rotation, latitude,
longitude, or frame transformation in M04.

| Quantity | Unit |
| --- | --- |
| Position | m |
| Velocity | m/s |
| Acceleration | m/s^2 |
| Time and timestep | s |
| Gravitational parameter \(\mu\) | m^3/s^2 |
| Specific orbital energy | m^2/s^2 = J/kg |
| Specific angular momentum | m^2/s |

## 4. How gravity and integration are connected

The M04 derivative does not repeat the gravity equation. It calls the canonical
M02 `two_body_acceleration()` function:

```text
Cartesian state
      |
      v
two-body state derivative
  dr/dt = velocity
  dv/dt = existing two_body_acceleration(position, mu)
      |
      v
existing Euler or RK4 single step
      |
      v
fixed-step propagation driver
      |
      v
time-stamped trajectory samples
```

The two-body model is autonomous: its derivative depends on state but not on
absolute time. The derivative API still accepts time because Euler/RK4 use the
general \(f(t,\mathbf x)\) interface; M04 deliberately leaves that argument
unused rather than pretending it affects gravity.

## 5. Fixed-step propagation policy

The propagation driver accepts a start time, end time, signed nominal timestep,
initial state, integration method, and derivative function. Its policy is:

- include the initial sample;
- require finite start/end times and a finite, nonzero timestep;
- require positive `dt` for forward propagation and negative `dt` for backward
  propagation;
- take full nominal steps while possible;
- shorten the final step to land exactly on the requested endpoint;
- never silently overshoot; and
- return one initial sample without evaluating the derivative when duration is
  zero.

The driver selects either `euler_step()` or `rk4_step()`. It does not contain a
second implementation of either integration formula.

## 6. Circular orbital speed and period

In a circular orbit, gravity supplies exactly the required centripetal
acceleration:

\[
\frac{v_c^2}{r}=\frac{\mu}{r^2}.
\]

Therefore

\[
v_c=\sqrt{\frac{\mu}{r}}.
\]

The orbital circumference is \(2\pi r\), so the period is

\[
T=\frac{2\pi r}{v_c}
=2\pi\sqrt{\frac{r^3}{\mu}}.
\]

Circular velocity is sideways—perpendicular to the radius. It is the speed at
which the spacecraft falls inward while moving sideways quickly enough for the
curved Earth to continually recede beneath it.

## 7. Why a satellite falls around Earth

Orbit does not mean there is no gravity. At 500 km altitude the point-mass
acceleration is about 8.43 m/s^2, still a large fraction of surface gravity.
Gravity constantly bends the velocity toward Earth. With insufficient sideways
speed the spacecraft intersects Earth; with the circular speed, its free-fall
path curves around Earth instead.

```text
strong inward gravitational acceleration
              +
large sideways velocity
              =
continuous free fall around Earth
```

## 8. Specific orbital energy

Specific mechanical energy is energy per unit spacecraft mass:

\[
\varepsilon
=
\frac{\lVert\mathbf v\rVert^2}{2}
-\frac{\mu}{\lVert\mathbf r\rVert}.
\]

Its units are m^2/s^2 or J/kg. For a circular orbit,

\[
\varepsilon=-\frac{\mu}{2r}.
\]

The ideal continuous two-body equations conserve this value. A numerical
trajectory whose energy changes is not solving those equations exactly.

## 9. Specific angular momentum

Specific angular momentum is

\[
\mathbf h=\mathbf r\times\mathbf v.
\]

Its units are m^2/s. The vector is normal to the orbital plane. Central gravity
creates no moment about Earth's center because
\(\mathbf r\times\mathbf a=\mathbf 0\), so both the direction and magnitude of
\(\mathbf h\) remain constant in the ideal continuous model.

## 10. Conservation laws as numerical diagnostics

An orbit can look plausible while accumulating serious numerical error.
Specific energy and angular momentum give quantitative, physics-based checks.
M04 records relative drift as

\[
\frac{\varepsilon(t)-\varepsilon_0}{|\varepsilon_0|}
\]

and

\[
\frac{|\mathbf h(t)|-|\mathbf h_0|}{|\mathbf h_0|}.
\]

The maximum absolute value over a trajectory summarizes worst observed drift.
These diagnostics measure how well the numerical trajectory respects the
selected ideal model; they do not establish that the model perfectly describes
the real Earth environment.

## 11. Worked 500 km circular orbit

The demo uses the same WGS 84 constants as M02:

```text
Earth reference radius       6,378,137 m
altitude                       500,000 m
orbital radius               6,878,137 m
Earth mu                3.986004418e14 m^3/s^2
```

The production helpers and demo calculate:

| Quantity | Value |
| --- | ---: |
| Gravitational acceleration | 8.425508710 m/s^2 |
| Circular speed | 7,612.608173 m/s = 7.612608 km/s |
| Analytical period | 5,676.978029 s = 94.616300 min |
| Specific orbital energy | -28,975,901.600 m^2/s^2 |
| Specific angular momentum magnitude | 52,360,561,942.754 m^2/s |

The initial state is

\[
\mathbf r_0=[6{,}878{,}137,0,0]\ \text{m},
\qquad
\mathbf v_0=[0,7{,}612.608173,0]\ \text{m/s}.
\]

Position and velocity are perpendicular, so this is an ideal equatorial
circular orbit in the simplified inertial frame. Inclination and orbital
elements are intentionally outside M04.

## 12. Why use a 10-second timestep?

Ten seconds is about 0.176% of the 94.6-minute orbital period. It produces
roughly 568 steps per orbit: inexpensive for this six-component state, fine
enough for excellent RK4 behavior, and still large enough to expose Forward
Euler's known long-duration weakness without choosing an absurd step merely to
make it fail.

## 13. Measured one-orbit behavior

Both methods use the same initial state, 10-second nominal timestep, and one
analytical period. The final interval is shortened so the endpoint is exactly
\(T\).

| Metric | Forward Euler | Classical RK4 |
| --- | ---: | ---: |
| Final position closure error | 4,152,340.133 m | 0.015806 m |
| Final velocity closure error | 4,086.812032 m/s | 1.74818e-5 m/s |
| Maximum radius deviation | 919,840.239 m | 0.003326 m |
| Maximum relative energy drift | 1.09568e-1 | 2.89584e-11 |
| Maximum relative angular-momentum drift | 5.96913e-2 | 1.44791e-11 |

Forward Euler evaluates only the beginning slope and advances position and
velocity with that stale information. In this coupled orbital problem it adds
energy and angular momentum, causing the numerical orbit to expand and fail to
close. RK4 samples how both velocity and acceleration change inside the step,
so its fourth-order truncation error is dramatically smaller.

## 14. Five-orbit accumulated drift

The same comparison over five analytical periods shows accumulated error:

| Metric | Forward Euler | Classical RK4 |
| --- | ---: | ---: |
| Final position closure error | 15,705,470.913 m | 0.097797 m |
| Maximum radius deviation | 3,312,113.815 m | 0.004123 m |
| Maximum relative energy drift | 3.14157e-1 | 1.44910e-10 |
| Maximum relative angular-momentum drift | 2.06992e-1 | 7.24550e-11 |

RK4 is much better here, but it is not exact. It still has finite-step
truncation and floating-point error, and classical RK4 is not a
conservation-preserving or symplectic method. Longer durations, different
steps, and different dynamics require fresh validation.

## 15. Physical-model error versus numerical error

These are distinct engineering error sources.

| Physical-model error | Numerical integration error |
| --- | --- |
| The chosen differential equation omits real physics. | The computer approximates the chosen equation. |
| M04 omits nonspherical gravity and J2, atmospheric drag, third-body gravity, solar radiation pressure, Earth rotation, and thrust. | M04 has Euler/RK4 truncation error, finite timestep effects, and floating-point roundoff. |
| Reducing `dt` does not restore omitted physics. | Reducing `dt` should reduce truncation error in the convergent regime. |
| A higher-fidelity force model addresses this source. | A suitable method, step policy, and numerical validation address this source. |

A trajectory can integrate the wrong physical model very accurately, or
integrate a useful model badly. Both questions must be evaluated separately.

## 16. Output and visualization

The C++ demo exports one-orbit and five-orbit CSVs containing:

```text
integrator, time_s,
x_m, y_m, z_m,
vx_mps, vy_mps, vz_mps,
radius_m, altitude_m, speed_mps,
specific_energy_m2_s2,
specific_angular_momentum_m2_s,
relative_energy_error,
relative_angular_momentum_error
```

Python reads these authoritative rows and performs only display conversions and
plotting. It contains no gravity or integration implementation. The generated
figures show the RK4 X-Y trajectory with equal axes, Euler/RK4 altitude, signed
relative energy error, and signed angular-momentum drift.

## 17. Assumptions, limitations, and failure modes

M04 assumes spherical point-mass Earth gravity, a negligible spacecraft mass,
an inertial Cartesian frame, and no perturbing forces. Important failure modes
include non-finite states, the gravity singularity at Earth's center, mixed
units or frames, a timestep with the wrong sign, an excessively large fixed
step, and interpreting visually plausible motion as sufficient validation.

There are no frame transformations, orbital elements, perturbations, attitude,
sensors, estimation, guidance, control, or rendezvous capabilities.

# What you should now understand

1. Why do we need both position and velocity to define translational orbital state?
2. Why is orbital motion naturally a second-order differential equation?
3. How do we rewrite it as a first-order system?
4. What does `dr/dt = v` mean physically?
5. What does `dv/dt = a` mean physically?
6. Why does a circular orbit require sideways velocity?
7. Why is a satellite still strongly affected by gravity in orbit?
8. What is specific orbital energy?
9. What is specific angular momentum?
10. Why should they remain constant in an ideal two-body simulation?
11. What does energy drift tell us about the numerical integrator?
12. Why does Euler perform badly for long orbital propagation?
13. Why is RK4 still not an exact solution?
14. What is the difference between numerical error and physical-model error?
15. Why might we eventually need adaptive integrators or higher-fidelity force models?
