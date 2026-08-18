# Lesson 005: Propagator Validation and Convergence

## Learning objective

M04 proved that AstraDock can propagate a circular orbit. M05 answers a more
important engineering question:

> **How accurate is the propagator, and how do we know when to trust it?**

```text
producing a trajectory that "looks right" ≠ trusting the numbers
measuring error ≠ trusting the numbers
characterizing convergence + invariants + regression = knowing when to trust
```

---

## 1. Why passing unit tests is not enough

M04's unit tests verify that the code implements the equations correctly.
But a simulation can implement the right equations with the wrong timestep,
and produce a trajectory that looks plausible while being numerically wrong.

Consider Forward Euler at $\Delta t = 10$ s for a 500 km orbit:
- The spacecraft "orbits" for the full period.
- The plot looks roughly circular.
- But the position error after one orbit is **4.1 million metres**.
- The energy has changed by **11%**.

A simulation that does not crash is not automatically trustworthy.

**Key insight:** Unit tests tell you whether the code is *correctly implemented*.
Validation tells you whether the *numerical solution is accurate enough* for
your purpose.

---

## 2. What is convergence?

A numerical method **converges** when reducing the timestep reduces the error.
If the error goes to zero as $\Delta t \to 0$, the method is convergent.

Not all methods converge. Not all methods converge at the same rate.

**Global truncation error** is the accumulated error at the end of the
integration. For a method of order $p$, the global error is roughly
proportional to $\Delta t^p$:

$$E \approx C \cdot \Delta t^p$$

where $C$ depends on the problem (derivative smoothness, integration interval).

---

## 3. Empirical convergence order

Given two timesteps $h_1 > h_2$ and their corresponding errors $E_1, E_2$,
the empirical order is:

$$p = \frac{\log(E_1 / E_2)}{\log(h_1 / h_2)}$$

For a half-step refinement ($h_2 = h_1 / 2$):
- If error halves $\to p \approx 1$ (first-order, like Euler)
- If error divides by 16 $\to p \approx 4$ (fourth-order, like RK4)

This is how you **measure** the actual order from your own simulation, rather
than relying on a textbook claim.

---

## 4. Measured M05 results

### RK4 on orbital dynamics

| $\Delta t$ (s) | Position error (m) | Empirical order |
| ---: | ---: | ---: |
| 40 | 4.77 | -- |
| 20 | 0.27 | 4.15 |
| 10 | 0.016 | 4.08 |
| 5 | 0.001 | 4.04 |
| 2.5 | 5.9e-05 | 4.02 |

RK4 is approximately **fourth-order** on the two-body orbital problem. As
$\Delta t$ shrinks, the measured order approaches exactly 4.

### Euler on orbital dynamics

| $\Delta t$ (s) | Position error (m) | Empirical order |
| ---: | ---: | ---: |
| 40 | 12,486,800 | -- |
| 20 | 7,502,956 | 0.73 |
| 10 | 4,152,340 | 0.85 |
| 5 | 2,187,391 | 0.92 |
| 2.5 | 1,122,635 | 0.96 |

Euler is approximately **first-order**. At the coarsest steps the measured
order is below 1 because the orbit spirals so far outward that the measured
"position closure error" is dominated by global trajectory deformation, not
just local truncation.

### Why the measured order is not exactly the theoretical order

The empirical order approaches the theoretical order as $\Delta t$ decreases.
At coarse steps, higher-order error terms contaminate the measurement. At very
fine steps, floating-point round-off ($\sim 10^{-16}$ for doubles) starts to
floor the error. The true asymptotic order lives in the middle.

---

## 5. Timestep sensitivity

Even a "correct" method can fail at a bad timestep:

- **Too coarse**: truncation error dominates, orbit spirals or fails to close.
- **Too fine**: floating-point round-off accumulates, CPU time wastes.
- **Just right**: truncation error is small but above the noise floor.

The M05 sweep shows that RK4 at $\Delta t = 10$ s is in the "just right" zone
for a 500 km orbit: 1.6 cm position error, $10^{-11}$ energy drift.

Euler at the same timestep is not: 4,152 km position error, 11% energy drift.

---

## 6. Truncation error vs floating-point error

| Error source | What it is | How it scales |
| --- | --- | --- |
| **Truncation** | The integrator approximates the ODE solution. Each step has a local error; errors accumulate globally. | $\sim \Delta t^p$ for order-$p$ method. Decreases as $\Delta t$ shrinks. |
| **Floating-point** | Each arithmetic operation loses $\sim 16$ decimal digits of precision. Errors accumulate over many operations. | $\sim \sqrt{N} \cdot \epsilon$ (random walk) or $\sim N \cdot \epsilon$ (systematic), where $N = T/\Delta t$. **Increases** as $\Delta t$ shrinks. |

The optimal timestep balances these two competing effects. For our 500 km
demonstration with RK4, truncation error dominates at all tested timesteps.

---

## 7. Phase error

A propagator can preserve the orbital radius while drifting in phase. This
means the spacecraft is at the right distance but the wrong angular position
along the orbit.

**How we measure it:**
1. Compute the unwrapped phase angle $\theta(t) = \operatorname{atan2}(y, x)$
   at each sample.
2. Compare with the analytical phase $\theta_0 + n t$, where $n = \sqrt{\mu/r^3}$.
3. Maximum absolute difference is the phase error.

Phase error matters for:
- **Rendezvous**: being at the wrong phase means being at the wrong place at
  the wrong time.
- **Ground track prediction**: phase drift means the satellite passes over a
  different longitude than expected.

RK4 at $\Delta t = 10$ s: phase error is 2.3 nanoradians — negligible.
Euler at $\Delta t = 10$ s: phase error is 0.56 radians — the spacecraft is
~32 degrees behind where it should be.

---

## 8. Energy drift as a diagnostic

In the ideal two-body model, specific orbital energy is conserved:
$$\epsilon = \frac{v^2}{2} - \frac{\mu}{r} = \text{constant}$$

A changing $\epsilon$ means the integrator is injecting or removing energy.
This is a **physics-based** diagnostic: it does not require knowing the
analytical trajectory. It only requires knowing what *should* be conserved.

| Integrator | $\Delta t$ | Max relative energy drift |
| --- | ---: | ---: |
| RK4 | 10 s | $2.9 \times 10^{-11}$ |
| Euler | 10 s | 0.11 |

Euler's 11% energy drift explains why the orbit spirals outward: it
consistently adds energy at each step.

---

## 9. Angular-momentum drift

Similarly, specific angular momentum $\mathbf{h} = \mathbf{r} \times \mathbf{v}$
should be constant in direction and magnitude for central gravity.

Magnitude drift tells you how well the integrator preserves the orbital
plane's "size." Direction drift tells you whether the orbital plane is
tilting numerically.

For the equatorial orbit used here, $\mathbf{h}$ points in the $+Z$ direction.
RK4 at $\Delta t = 10$ s keeps the direction drift below $10^{-8}$ radians.

---

## 10. Analytical references

The best validation compares against a **known** answer. For a circular orbit:

| Analytical quantity | Formula |
| --- | --- |
| Circular speed | $v_c = \sqrt{\mu/r}$ |
| Period | $T = 2\pi\sqrt{r^3/\mu}$ |
| Specific energy | $\epsilon = -\mu/(2r)$ |
| Angular momentum | $h = r \cdot v_c$ |
| Mean motion | $n = \sqrt{\mu/r^3}$ |

These give exact reference values against which to compare numerical results.
**They are exact** for the chosen idealized model.

---

## 11. Numerical reference solutions

Sometimes no analytical solution exists. In that case, a **numerical reference**
is produced by running the same integrator at a much finer timestep:

- RK4 at $\Delta t = 0.25$ s $\approx$ "high-accuracy numerical reference"
- Not the same as exact truth — still has truncation and floating-point error
- Document it explicitly as a numerical reference, not mathematical truth

---

## 12. Regression tests

A regression test records a known-good numerical result and checks that future
code changes do not degrade it:

- RK4, $\Delta t = 10$ s, one orbit, position closure $< 0.1$ m
- RK4, $\Delta t = 10$ s, one orbit, energy drift $< 10^{-9}$
- Euler, $\Delta t = 10$ s, one orbit, position error $> 100$ km (sanity check)

These are not "unit tests" in the pure sense — they propagate a full orbit.
But they are deterministic and fast enough to run on every build.

---

## 13. Verification vs Validation

These are distinct engineering concepts:

**Verification:** "Did we implement the equations/algorithm correctly?"

- Unit tests: does `two_body_acceleration()` return the right value?
- Analytical comparison: does RK4 match the known analytical solution?
- Invariant preservation: is energy conserved (within truncation)?
- Code reviews: is the implementation free of bugs?

**Validation:** "Does the model and simulation behave appropriately for the
problem we are trying to represent?"

- Does the two-body model capture the physics we need?
- Is the timestep small enough for our accuracy requirements?
- Would a real spacecraft behave this way?

**For AstraDock at M05:**
- Unit tests = verification
- Analytical comparison = verification
- Invariant preservation = verification / diagnostic evidence
- Comparison against a higher-fidelity physical reference = future validation

---

## 14. Why RK4 is not symplectic

A **symplectic** integrator preserves the geometric structure of Hamiltonian
mechanics. For orbital mechanics, this means:
- Energy oscillates around the true value but does not drift monotonically.
- Long-duration orbits stay stable even with moderate timesteps.

RK4 is **not** symplectic. Over very long durations (thousands of orbits),
it still accumulates a small but systematic energy drift. Symplectic methods
(Verlet, Gauss-Legendre) are preferred for very long integrations, but RK4 is
excellent for the moderate-duration scenarios we currently need.

---

## 15. Why an accurate short simulation can still accumulate long-term error

M04 showed that RK4 at $\Delta t = 10$ s has 1.6 cm error after one orbit.
Over 1,000 orbits, the error might grow to meters or more. Error accumulation
is not always linear — it depends on the dynamics. But the principle is:

**Validating for one orbit does not prove accuracy for a thousand orbits.**

---

## 16. Why physical-model error differs from numerical error

| Aspect | Numerical error | Physical-model error |
| --- | --- | --- |
| Source | Approximating the ODE solution | The ODE omits real physics |
| Fix | Smaller timestep, better integrator | Better force model (J2, drag, ...) |
| Measurement | Compare timesteps, check invariants | Compare against higher-fidelity model |
| M05 addresses | Yes | No — by design |

You can integrate the wrong model perfectly accurately. That makes a very
precise but physically wrong trajectory.

---

# What you should now understand

1. **Why can a simulation produce a visually plausible orbit while still being
   numerically inaccurate?** The orbit looks roughly circular even when the
   position is millions of metres wrong. Visual inspection is not validation.

2. **What does convergence mean?** As $\Delta t \to 0$, the numerical solution
   approaches the true solution of the ODE. The rate of convergence is the
   method's order $p$.

3. **How can we estimate numerical order from timestep refinement?**
   $p = \log(E_1/E_2) / \log(h_1/h_2)$. Measure errors at two timesteps and
   compute the ratio of log-errors.

4. **Why should RK4 approach fourth-order convergence rather than equal
   exactly 4 at every timestep?** At coarse steps, higher-order error terms
   contaminate the measurement. At very fine steps, floating-point noise
   floors the error. The asymptotic order $p = 4$ lives in between.

5. **Why does Euler inject large orbital energy error?** It uses the slope at
   the beginning of each step for the entire interval. In orbital dynamics,
   both position and velocity are coupled — using a stale slope means
   consistently overestimating or underestimating the true trajectory, which
   manifests as energy injection.

6. **What is phase error?** The difference between where the numerical solution
   is along the orbit and where the analytical solution says it should be. It
   measures angular lag/lead, not radial error.

7. **Why can small radial error coexist with significant phase error?** An
   integrator can keep the spacecraft at approximately the right distance from
   Earth while running slightly fast or slow along the orbit. The radius is
   nearly constant, but the spacecraft arrives at the wrong point on time.

8. **Why are energy and angular momentum useful diagnostics?** They are
   conserved quantities in the ideal two-body model. A changing value means
   the integrator is not solving the equations exactly. They provide
   physics-based checks that do not require a reference trajectory.

9. **What is the difference between verification and validation?** Verification
   checks whether the code implements the equations correctly. Validation
   checks whether the model and numerical solution are appropriate for the
   problem being studied.

10. **Why is a very small timestep not automatically the best engineering
    choice?** As $\Delta t$ shrinks, truncation error decreases but floating-point
    round-off and computation time increase. The optimal timestep balances
    these effects for the required accuracy.

11. **Why is RK4 not necessarily ideal for arbitrarily long orbital
    propagation?** RK4 is not symplectic — it does not preserve the geometric
    structure of Hamiltonian systems. Over many thousands of orbits, it can
    accumulate systematic energy drift. Symplectic integrators are preferred
    for very long integrations.

12. **What would we need to model to make Earth dynamics more realistic?**
    J2 (Earth's oblateness), atmospheric drag, third-body gravity (Moon, Sun),
    solar radiation pressure, Earth rotation, and possibly general-relativistic
    corrections for high-precision applications.