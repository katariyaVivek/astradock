# Lesson 003: Numerical Integration

## The engineering problem

A differential equation describes how a quantity changes rather than giving
its future value directly. A first-order ordinary differential equation (ODE)
has the form

\[
\frac{d\mathbf{x}}{dt}=f(t,\mathbf{x}).
\]

Here:

- \(t\) is the current time,
- \(\mathbf{x}\) is the current **state**,
- \(f(t,\mathbf{x})\) is the **derivative function**, and
- \(\Delta t\) is the finite time step over which we want to advance.

The state contains the information needed to continue the problem. It may be a
single number, such as temperature, or a vector of values. The derivative has
state-units per unit time and reports how quickly that state is changing at one
instant.

A simulator usually knows

```text
current state + how quickly the state is changing
```

but does not know the future state directly. Numerical integration estimates
that future state by sampling the derivative over a finite interval.

For spacecraft translation, a later milestone will use the conceptual
relationships

\[
\frac{d\mathbf r}{dt}=\mathbf v,
\qquad
\frac{d\mathbf v}{dt}=\mathbf a.
\]

M03 does not assemble that position-and-velocity state and does not call the
two-body gravity model. The curriculum boundary is deliberate:

```text
M02: physics supplies an instantaneous derivative
M03: an integrator advances any suitable state from a derivative
M04: combine those ideas for orbital propagation
```

## Frames and units

Euler and RK4 are mathematical algorithms, so they do not assign a coordinate
frame or physical unit to the state. Those meanings belong to the caller and
the state definition. If \(t\) and \(\Delta t\) are seconds, then the derivative
must be expressed in state-units per second. A `Vector3` used as a state must
still have its frame and component units documented at the surrounding API.

This separation is useful: the same integration code can advance a scalar
decay problem, a frame-documented vector, or a future structured spacecraft
state without pretending those states mean the same thing.

## Forward Euler

Forward Euler uses the derivative at the beginning of the interval:

\[
\mathbf{x}_{n+1}
=
\mathbf{x}_n+\Delta t\,f(t_n,\mathbf{x}_n).
\]

Its assumption is simple:

> Treat the current derivative as constant during the entire step.

Geometrically, Euler follows the tangent line from the current solution point.
It requires only one derivative evaluation, is easy to inspect, and exposes
the central idea of numerical integration clearly. It is therefore valuable
for learning and as a comparison baseline. The same simplicity also means it
can miss curvature in the true solution and accumulate substantial error.

### A manual Euler step

Consider

\[
y'=-y,\qquad y(0)=1,\qquad \Delta t=0.1.
\]

At the beginning, the derivative is \(y'(0)=-1\). Euler gives

\[
y_1=1+0.1(-1)=0.9.
\]

The analytical value is

\[
y(0.1)=e^{-0.1}\approx 0.904837418.
\]

The one-step absolute error is therefore about \(0.004837418\). The exact
solution curves during the interval because its slope becomes less negative;
Euler used only the initial, more-negative slope.

## Classical fourth-order Runge--Kutta

Classical RK4 samples four slopes:

\[
k_1=f(t_n,\mathbf{x}_n),
\]

\[
k_2=f\left(t_n+\frac{\Delta t}{2},
\mathbf{x}_n+\frac{\Delta t}{2}k_1\right),
\]

\[
k_3=f\left(t_n+\frac{\Delta t}{2},
\mathbf{x}_n+\frac{\Delta t}{2}k_2\right),
\]

\[
k_4=f(t_n+\Delta t,\mathbf{x}_n+\Delta t\,k_3),
\]

then combines them:

\[
\mathbf{x}_{n+1}
=
\mathbf{x}_n
+\frac{\Delta t}{6}(k_1+2k_2+2k_3+k_4).
\]

The stages estimate the slope at the start, at two independently predicted
midpoints, and at a predicted endpoint. The weighted average accounts for how
the derivative changes within the step. RK4 therefore uses four derivative
evaluations where Euler uses one, trading additional computation for much
better accuracy on suitable smooth problems.

## Local error, global error, and method order

**Local truncation error** is the error introduced by one numerical step when
that step begins from the exact state. **Global error** is the accumulated
error at a fixed final time after taking many steps; each new step begins from
an already approximate state.

For a sufficiently smooth ODE:

| Method | Local truncation error | Global error over a fixed interval |
| --- | --- | --- |
| Forward Euler | \(O(\Delta t^2)\) | \(O(\Delta t)\) |
| Classical RK4 | \(O(\Delta t^5)\) | \(O(\Delta t^4)\) |

Euler is therefore a **first-order** method. Once the step is small enough for
the leading error term to dominate, halving \(\Delta t\) should approximately
halve its global error. RK4 is a **fourth-order** method, so halving the step
should reduce its global error by approximately \(2^4=16\). These are
asymptotic convergence expectations, not exact identities for every equation,
step size, or floating-point environment.

## Error-versus-step-size experiment

The M03 tests integrate

\[
y'=-y,\qquad y(0)=1
\]

from \(t=0\) to \(t=1\), then compare with \(e^{-1}\). The following values are
computed from the same fixed-step recurrences exercised by the tests:

| \(\Delta t\) | Euler absolute error | RK4 absolute error |
| ---: | ---: | ---: |
| 0.2 | 4.019944117e-2 | 5.796953860e-6 |
| 0.1 | 1.920100107e-2 | 3.332410564e-7 |
| 0.05 | 9.393518763e-3 | 1.997609661e-8 |
| 0.025 | 4.647001284e-3 | 1.222742074e-9 |

For successive halvings, the Euler error ratios are approximately 2.09, 2.04,
and 2.02. The RK4 ratios are approximately 17.40, 16.68, and 16.34. Both trends
approach their theoretical values. At \(\Delta t=0.1\), RK4's final error is
over fifty thousand times smaller for this particular smooth problem. That is
evidence for this problem and step size, not a universal performance promise.

The tests also integrate the vector ODE

\[
\frac{d\mathbf{x}}{dt}=-\mathbf{x},
\qquad
\mathbf{x}(0)=[1,2,-3],
\]

whose analytical solution is

\[
\mathbf{x}(t)=e^{-t}[1,2,-3].
\]

This demonstrates that the implementation applies the same mathematics to
`Vector3` and is not accidentally scalar-only.

## Convergence, cost, and accumulated numerical error

Reducing \(\Delta t\) generally improves an accurate, stable discretization
because each straight or polynomial approximation spans a shorter interval.
It also increases cost: halving the step doubles the number of steps over a
fixed duration. Each Euler step evaluates the derivative once, while each RK4
step evaluates it four times.

Smaller steps also mean more floating-point operations. At first truncation
error normally falls, but at extremely small steps roundoff and cancellation
can become relevant. Real simulations choose a step by measuring both accuracy
and cost against requirements rather than assuming "smaller" is sufficient.

## Stability is separate from order

Convergence order describes how error changes as the step approaches zero; it
does not guarantee that a chosen finite step behaves well. For the test problem
\(y'=-y\), Euler advances by the factor \(1-\Delta t\). If a positive step is
too large, the numerical solution can oscillate or grow even though the true
solution decays.

RK4 has a larger useful stability region than Euler for many problems, but it
is still an explicit fixed-step method with limits. A high-order method can be
unstable on a stiff problem or at an unsuitable step size. Aerospace
simulators must therefore validate accuracy, convergence, stability, and
conserved or bounded quantities for the actual dynamics and operating range.

## API and time-step policy

M03 provides two header-only function templates:

```cpp
template <typename State, typename Derivative>
State euler_step(double t, const State& state, double dt, Derivative&& derivative);

template <typename State, typename Derivative>
State rk4_step(double t, const State& state, double dt, Derivative&& derivative);
```

The state needs only the addition and scalar-multiplication operations used by
the selected algorithm. This small amount of genericity is justified because
both `double` and `Vector3` are verified now, and M04 can introduce a
translational state without duplicating Euler and RK4. There is no integrator
base class, inheritance hierarchy, plugin system, or adaptive-step framework.

The deterministic time policy is:

- `t` and `dt` must be finite;
- the midpoint and endpoint times used by the step must remain representable;
- `dt == 0` returns the unchanged state without evaluating the derivative;
- negative `dt` is supported for backward integration; and
- positive `dt` is the normal forward-in-time case.

For the M03 state types, `double` and `Vector3`, non-finite input states,
derivatives, intermediate states, and results are rejected. Invalid supplied
values and derivative outputs raise `std::domain_error`; arithmetic that
produces a non-representable result raises `std::overflow_error`. Exceptions
raised by the derivative function itself propagate to the caller.

Responsibility remains layered:

- the **integrator** implements the step equations and validates time plus the
  known M03 numeric states;
- the **state type** supplies meaningful arithmetic and must define validation
  for any custom internal components;
- the **derivative function** supplies the correct mathematical model and
  compatible units; and
- the future **simulation layer** will choose step size, duration, sampling,
  frames, error budgets, and scenario-level validation.

## Known limitations

These primitives are fixed-step, explicit, and single-step. They provide no
adaptive error estimate, dense output, event detection, stiffness detection,
automatic step selection, or propagation loop. M03 has intentionally not
created a position/velocity spacecraft state and has not connected the
integrators to `two_body_acceleration()`.

# What you should now understand

1. What does an ODE describe?
2. What information does the derivative function return?
3. What assumption does Forward Euler make during each step?
4. Why does decreasing `dt` usually reduce Euler error?
5. What does first-order convergence mean?
6. Why does RK4 evaluate four slopes?
7. What does fourth-order convergence mean?
8. Why is RK4 usually more accurate than Euler at the same `dt`?
9. Why does higher order not automatically guarantee stability?
10. Why should we validate the integrator before connecting it to orbital dynamics?
