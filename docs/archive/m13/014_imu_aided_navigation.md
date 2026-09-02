# Lesson 014 — IMU-Aided Translational Navigation

## 1. The Physical Problem: Propagating Without a Model Crutch

M13A predicted between GNSS fixes by re-integrating its own dynamics model
(two-body gravity). That is legitimate only when reality actually obeys the
model. A real spacecraft cannot assume that: atmospheric drag, thruster
firings, and gravity-model errors all make the true acceleration differ from
the modeled one. The spacecraft must therefore **sense** how it accelerates.
That is the job of an Inertial Measurement Unit (IMU):

```text
accelerometer: measures specific force  f  (BODY frame)
gyroscope:     measures angular velocity w (BODY frame)
```

M13B replaces dynamics-only prediction with **inertial dead reckoning**:

```text
IMU measurement -> prediction -> GNSS update -> posterior estimate
        (100 Hz)                     (1 Hz)
```

## 2. Specific Force Is Not Acceleration

The single most misunderstood sensor in aerospace:

```text
f = a - g        (specific force = inertial acceleration minus gravity)
```

An accelerometer is a force-referenced instrument: internally a proof mass on
springs. Gravity acts on the proof mass AND the case equally, so springs stay
relaxed under free fall. Consequences:

- Sitting on the launch pad, f = +1 g upward (the pad pushes you).
- In orbit, g is large (~8 m/s²) but f ≈ **0** — the spacecraft is falling.
- Drag or thrust produce tiny nonzero readings (drag at 500 km ~ 1 µm/s²).

So to recover the acceleration that drives translation, the filter must undo
the measurement physics — this is the central equation chain of M13B:

```text
accelerometer reading f_m (BODY)
      |  subtract known bias   b_a
      v
corrected specific force (BODY)
      |  rotate with ATTITUDE ESTIMATE q_hat
      v
specific force (ECI)
      |  add gravity back       g(r_hat)
      v
inertial acceleration a_I (ECI)  -> integrate -> position + velocity
```

Production implementation (one function, independently audited):

```cpp
a_I = C_I_B(q_hat) * (f_m - b_a) + g(r_hat);
```

**Why gravity must be added back**: the accelerometer *removed* it (f = a − g),
so integration only makes sense after restoring it. In orbit this step looks
trivial (f ≈ 0 so a ≈ g) but it is conceptually everything: get the sign or
frame wrong here and the filter integrates the wrong world.

## 3. Attitude Is the Bridge Between Frames

`C_I_B(q_hat)` rotates body coordinates into ECI coordinates. The rotation is
only as good as the attitude estimate:

```text
delta-attitude
    -> wrong BODY->ECI rotation
    -> wrong inertial acceleration (error ~ |f| * delta-theta)
    -> velocity error
    -> position error
```

Measured in the M13B diagnostic test: a 10 mrad attitude error applied to a
0.1 m/s² specific force corrupts the reconstructed acceleration by
2·sin(δθ/2)·|f| ≈ 1 mm/s² — **larger than a typical MEMS accelerometer noise
floor**. In free fall |f| ≈ 0, so translational navigation tolerates poor
attitude during orbital coast; under sustained thrust or drag it does not.

> **M13B assumes the attitude is externally supplied/known.** The filter
> receives it through an explicit navigation-state input
> (`NavigationAttitudeEstimate`) — never as truth access. This is deliberate
> staging: M13C will estimate attitude jointly and fold its uncertainty into
> the covariance.

## 4. Multi-Rate Prediction/Update

Two sensors, two rates, one filter:

```text
t_k (every 10 ms):  IMU sample  -> predict_with_imu(dt = t_k - t_{k-1})
t_j (every 1 s):    GNSS fix    -> update_gnss(z)
```

Rules learned here:

- **Timestamp discipline**: `dt` comes from actual measurement timestamps,
  never an assumed constant. Binary floating point makes even "exact" 10 ms
  schedules wobble by nanoseconds; the filter honors them.
- **Zero-order hold (ZOH)**: within one prediction interval both the measured
  specific force and the supplied attitude are held constant; gravity varies
  continuously through RK4. At 100 Hz the ZOH error is far below sensor noise.
- Between GNSS fixes the estimate is pure dead reckoning; the fix snaps the
  estimate back and shrinks the covariance.

## 5. Process Noise Comes From the Sensor, Not From Tuning

M13A's Q was driven by an abstract unmodeled-acceleration sigma. In M13B the
dominant stochastic forcing between updates IS the accelerometer noise, so Q
is derived directly from it using the exact continuous-white-noise
discretization of the double integrator:

```text
Q = sigma_f^2 [ dt^3/3 I   dt^2/2 I ]
              [ dt^2/2 I   dt     I ]
```

Chained over N steps this reproduces the continuous result
Var[position] = sigma_f² T³/3 — the discretization is consistent, not ad hoc.
No "turn Q until the filter looks good" allowed: sigma_f is a physical
property of the instrument.

## 6. Why an IMU Cannot Tell You Position Forever

Dead reckoning chains integrations:

```text
accelerometer -> (integrate) -> velocity -> (integrate) -> position
```

Integration accumulates every constant error:

```text
tiny bias b
   -> velocity error grows LINEARLY:   delta_v = b*t
   -> position error grows QUADRATICALLY: delta_r = b*t^2 / 2
```

The M13B bias experiment (identity attitude, b along one axis, t = 300 s):

| true bias | position drift | velocity drift | 0.5·b·t² |
|---|---|---|---|
| 0 | ~0 | ~0 | 0 |
| 1e-4 m/s² | 4.58 m | 0.031 m/s | 4.50 m |
| 5e-4 m/s² | 22.89 m | 0.155 m/s | 22.50 m |

Drift is exactly linear in bias and quadratic in time (the small excess above
the flat-space value is radial gravity-gradient amplification). This table is
the entire economic argument for bias states in M13C/D: an unestimated bias
that is invisible on a 10-second horizon destroys navigation on a 30-minute
horizon. GNSS fixes hide the drift while they arrive — the dropout experiment
shows what happens when they stop.

## 7. Dropout: Dead Reckoning Under Pressure

With GNSS unavailable for 100 s the filter runs prediction-only:

```text
position sigma : 1.22 m -> 2.26 m  (monotone growth, no measurements)
max position error during outage: 1.56 m   (dead reckoning holds the line)
after recovery  : sigma contracts to 1.29 m within one update cycle
```

Covariance growth during information starvation and contraction when
measurements return are THE qualitative signatures of a healthy filter. Not
every component shrinks monotonically after recovery — geometry and process
noise shape the trajectory of P; only gross monotone behavior is guaranteed.

## 8. What the Gyroscope Is For (and Isn't Yet)

The gyro sample flows through the same pipeline and timestamps as the
accelerometer, but M13B does NOT use it inside the filter: attitude
estimation is deferred to M13C. Building an ad-hoc attitude integrator inside
M13B would duplicate M09/M10 machinery behind the estimator's back and create
a second attitude-estimation path — forbidden. Honest staging beats fake
completeness.

## 9. Verification Culture

Three layers, all required before believing the implementation:

1. **Analytical**: hand-computed DCM direction checks (90° about z maps body
   x̂ to ECI ŷ — a conjugated DCM maps it to −ŷ), independent textbook-DCM
   reconstruction across nontrivial attitudes, first-order small-dt limits,
   hand-built covariance propagation, finite-difference audits of the
   dynamics gradient across attitudes.
2. **Simulation/statistical**: bitwise degeneracy to the M13A predictor when
   f − b = 0, free-fall tracking against the open-loop reference at nm scale,
   drag capture through the sensor (76 m coast divergence vs sub-metre aided),
   dropout growth/recovery, determinism, truth non-interference, NEES/NIS
   Monte Carlo consistency.
3. **Independent oracle**: a from-scratch Python IMU-EKF replaying recorded
   measurements agrees with C++ to 4.4e-16 m — machine precision between two
   independent implementations of the same mathematics.

Debugging war story worth remembering: the drag scenario initially showed the
IMU-aided filter performing no better than coasting. Root cause was NOT the
filter — the harness measured specific force in the frame of one attitude
trajectory while rotating it back with another. Same-frame consistency
between the sensor model and the navigation input is part of the known-
attitude contract.

## What You Should Now Understand

1. Why does an accelerometer read ZERO in orbit but 1 g on the pad?
2. Write the central reconstruction equation and name each term's frame.
3. Where does attitude enter, and how does attitude error scale into
   position error? Why is free fall forgiving?
4. Why must Q derive from the accelerometer noise instead of being tuned?
5. Derive (or justify) the dt³/dt²/dt structure of the white-noise-acceleration Q.
6. How does a constant accelerometer bias propagate into velocity and
   position error, and why does quadratic growth make long missions hopeless
   without estimation?
7. What happens to the estimate and covariance during a GNSS outage?
8. Why can't M13B use the gyro inside the filter, and what will M13C change?
9. Why does ZOH on specific force suffice at 100 Hz?
10. What three verification layers does M13B require before trusting results?
