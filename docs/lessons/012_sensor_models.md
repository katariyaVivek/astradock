# Lesson 012 — Spacecraft Sensor Simulation & Measurement Models

## 1. The Physical Problem: Truth vs. Measurement

In classical flight simulations (Milestones M01–M11), we have direct access to the **truth state** of the spacecraft:
$$\mathbf{x}_{\text{truth}}(t) = \begin{bmatrix} \mathbf{r}_{\text{ECI}}(t) \\ \mathbf{v}_{\text{ECI}}(t) \\ q_{\mathcal{I}\_\mathcal{B}}(t) \\ \boldsymbol{\omega}_{\mathcal{B}}(t) \end{bmatrix}$$

In real aerospace missions, the onboard flight computer **never has direct access to the truth state**. Instead, the navigation system must infer the spacecraft state from a collection of imperfect, noisy, biased, rate-limited, and occasionally failing **sensor measurements**:

```text
       Truth Spacecraft Dynamics (ECI / BODY)
                         │
                         ▼
             Physical Sensor Hardware
    (Optical cameras, gyros, quartz proof masses, RF)
                         │
    ┌────────────────────┼────────────────────┐
    ▼                    ▼                    ▼
Imperfections        Cadence             Failures
• Constant Biases    • High (IMU 100Hz)   • Solar blinding
• White Noise        • Med (Star 10Hz)    • Earth occultation
• Misalignment       • Low (GNSS 1Hz)     • Signal dropouts
    │                    │                    │
    └────────────────────┼────────────────────┘
                         ▼
        Discrete Sensor Measurement Packets
    (Timestamp, Valid Flag, Measured Values)
                         │
                         ▼
       Future Navigation Filter / EKF (M13+)
```

---

## 2. Critical Concept: Why an Accelerometer Reads Nearly Zero in Orbit

One of the most common misconceptions in aerospace dynamics is the belief that an accelerometer measures gravitational acceleration $\mathbf{g}$.

### The Governing Principle of Specific Force
An accelerometer measures **specific force** $\mathbf{f}$—the net non-gravitational contact force per unit mass exerted on the sensor proof mass by its suspension:
$$\mathbf{f} = \mathbf{a}_{\text{inertial}} - \mathbf{g}(\mathbf{r})$$

where:
- $\mathbf{a}_{\text{inertial}} = \ddot{\mathbf{r}}$ is the total kinematic acceleration of the spacecraft relative to an inertial reference frame (ECI).
- $\mathbf{g}(\mathbf{r}) = -\frac{\mu}{r^3}\mathbf{r} + \mathbf{a}_{J2} + \dots$ is the local gravitational acceleration vector.

### Orbital Free Fall
In an ideal unperturbed orbit governed solely by gravity, Newton's second law gives:
$$\mathbf{a}_{\text{inertial}} = \mathbf{g}(\mathbf{r})$$

Substituting this into the accelerometer definition:
$$\mathbf{f}_{\text{ideal}} = \mathbf{g}(\mathbf{r}) - \mathbf{g}(\mathbf{r}) \equiv \mathbf{0}$$

> **Key Takeaway**: A spacecraft in circular orbit at $500\,\text{km}$ altitude experiences a strong gravitational acceleration of $g \approx 8.43\,\text{m/s}^2$. However, because the entire spacecraft and its proof mass are in continuous free fall together, the **ideal specific force is identically zero** ($f \approx 0$).

### Non-Gravitational Forces Produce Specific Force
When non-gravitational forces act on the spacecraft bus—such as atmospheric drag $\mathbf{a}_{\text{drag}}$, solar radiation pressure, or thruster firings $\mathbf{F}_{\text{thrust}}/m$—these contact forces accelerate the spacecraft relative to the free-falling frame:
$$\mathbf{f}_I = \mathbf{a}_{\text{non\_grav}, I} + \frac{\mathbf{F}_{\text{thrust}, I}}{m}$$

Transforming this inertial specific force into the spacecraft **Body frame** ($\mathcal{B}$) via the attitude quaternion:
$$\mathbf{f}_B = q^* \otimes \mathbf{f}_I \otimes q \equiv C_{\mathcal{B}\_\mathcal{I}} \mathbf{f}_I$$

The accelerometer measurement model is therefore:
$$\mathbf{f}_{\text{meas}} = \mathbf{f}_B + \mathbf{b}_a + \mathbf{n}_a \quad [\text{m/s}^2, \text{ BODY}]$$

---

## 3. Sensor Suite Mathematical Models

AstraDock models four core spacecraft sensors, each with explicit frame definitions, SI units, update cadences, and error statistics:

| Sensor | Physical Quantity | Measurement Frame | Units | Update Rate | Error Model |
| :--- | :--- | :--- | :---: | :---: | :--- |
| **IMU Gyroscope** | Body angular rate $\boldsymbol{\omega}_B$ | Spacecraft **BODY** ($\mathcal{B}$) | $\text{rad/s}$ | $100\,\text{Hz}$ | $\boldsymbol{\omega}_{\text{meas}} = \boldsymbol{\omega}_{\text{true}} + \mathbf{b}_g + \mathbf{n}_g$ |
| **IMU Accelerometer** | Specific force $\mathbf{f}_B$ | Spacecraft **BODY** ($\mathcal{B}$) | $\text{m/s}^2$ | $100\,\text{Hz}$ | $\mathbf{f}_{\text{meas}} = \mathbf{f}_B + \mathbf{b}_a + \mathbf{n}_a$ |
| **GNSS Receiver** | Absolute position & velocity | Inertial **ECI** ($\mathcal{I}$) | $\text{m}$, $\text{m/s}$ | $1\,\text{Hz}$ | $\mathbf{r}_{\text{meas}} = \mathbf{r}_{\text{true}} + \mathbf{b}_r + \mathbf{n}_r$, $\mathbf{v}_{\text{meas}} = \mathbf{v}_{\text{true}} + \mathbf{b}_v + \mathbf{n}_v$ |
| **Star Tracker** | Attitude quaternion $q_{\mathcal{I}\_\mathcal{B}}$ | **ECI from BODY** | dimensionless | $10\,\text{Hz}$ | $q_{\text{meas}} = (\delta q \otimes q_{\text{bias}} \otimes q_{\text{true}}).normalized()$ |
| **Range Sensor** | Relative range $\rho$ | Scalar Euclidean norm | $\text{m}$ | $10\,\text{Hz}$ | $\rho_{\text{meas}} = \max(0, \|\mathbf{r}_{\text{target}} - \mathbf{r}_{\text{sc}}\| + b_\rho + n_\rho)$ |

---

## 4. Star Tracker: Proper SO(3) Attitude Perturbations

Quaternions represent rotations on the 3-sphere $S^3$. A common novice error is adding unconstrained Gaussian noise directly to the quaternion four-vector ($q_{\text{bad}} = q_{\text{true}} + \mathbf{n}$). This violates the physical geometry of $SO(3)$ and destroys the group structure.

In AstraDock, attitude perturbations are generated strictly as **physical 3D rotations**:
1. Draw a random unit rotation axis $\hat{\mathbf{u}} \in S^2$ uniformly distributed on the unit sphere.
2. Draw an angular perturbation angle $\theta_{\text{err}} \sim \mathcal{N}(0, \sigma^2)$ from a zero-mean normal distribution with 1-sigma standard deviation $\sigma$ (e.g. $0.2\,\text{mrad} \approx 41\,\text{arcsec}$).
3. Form the perturbation quaternion $\delta q = \left[ \cos(\theta_{\text{err}}/2), \; \hat{\mathbf{u}} \sin(\theta_{\text{err}}/2) \right]$.
4. Apply the rotation perturbation via Hamilton quaternion multiplication:
   $$q_{\text{meas}} = (\delta q \otimes q_{\text{bias}} \otimes q_{\text{true}}).normalized()$$

This guarantees:
- Exact unit norm preservation ($\|q_{\text{meas}}\| = 1$).
- Double-cover invariance ($q$ and $-q$ represent the same physical attitude).
- Meaningful geodesic angular error metrics ($\theta_{\text{err}} = 2\arccos(|q_{\text{meas}} \cdot q_{\text{true}}|)$).

---

## 5. Why Sensor Update Rates Matter (Multi-Rate Navigation)

Spacecraft sensors operate at vastly different physical bandwidths and technological frequencies:

```text
Time (s) ───┼───────┼───────┼───────┼───────┼───────┼───────┼───────┼───────┼───►
IMU (100Hz) |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Star (10Hz) ●       ●       ●       ●       ●       ●       ●       ●       ●
Range (10Hz)▲       ▲       ▲       ▲       ▲       ▲       ▲       ▲       ▲
GNSS (1Hz)  ■                                                       ■
```

1. **IMU (100 Hz)**: Provides rapid, continuous kinematic integration between orbital and attitude updates. However, uncorrected IMU bias causes rapid quadratic error drift over time.
2. **Star Tracker (10 Hz)**: Provides high-accuracy absolute attitude reference to bound gyroscope bias drift.
3. **GNSS (1 Hz)**: Provides absolute position and velocity fixes to periodically correct translational orbital dead-reckoning.

In **Milestone M13**, the Extended Kalman Filter (EKF) will utilize this exact multi-rate structure: propagating state estimates at high frequency using the IMU ($100\,\text{Hz}$), and executing discrete Kalman measurement updates whenever lower-rate sensors ($10\,\text{Hz}$ Star Tracker, $1\,\text{Hz}$ GNSS) produce a valid sample.

---

## 6. Deterministic Stochastic Simulation & Dropouts

### Deterministic Random Number Generation
All stochastic models in AstraDock use an explicit, self-contained `DeterministicRng` wrapping 64-bit Mersenne Twister (`std::mt19937_64`). No global random singletons are permitted. Given the same random seed, simulations produce identical telemetry on every run.

### Sensor Dropouts & Validity Flags
Sensors can experience temporary signal loss (e.g. star tracker blinded by Earth limb or Sun, GNSS antenna occlusion). AstraDock models this through deterministic `DropoutWindow` intervals:
```cpp
struct DropoutWindow {
    double start_time_s;
    double end_time_s;
};
```
During an active dropout window:
- The measurement packet is marked with `valid = false`.
- Measured vector components are cleared to zero.
- The truth dynamics simulation continues completely undisturbed (**Truth Non-Interference Contract**).
