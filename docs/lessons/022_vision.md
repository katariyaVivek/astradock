# Lesson 022 — Computer Vision & Optical Navigation

## 1. Physical Problem & Motivation

GNSS gives meters; docking needs centimeters. A docking camera watching known
fiducials on the target closes that gap optically: 3-D body points → 2-D
pixels through pinhole optics, then pixels → relative pose through geometric
optimization (PnP). This lesson is classical vision FIRST (M23 learns only
against this baseline): no networks, no features-learned-end-to-end — just
projective geometry, least squares, and honest failure reports.

## 2. Camera Model (Four Frames)

WORLD (target body, m) → CAMERA (+Z boresight, +X right, +Y down) via
`p_c = R(p_w − t)` → NORMALIZED (`x_n = X/Z`) → PIXEL (`u = fx·x_d + cx` with
radial distortion `s = 1 + k1·r²`). Behind-clip (`Z ≤ 0.1 m`) and
out-of-bounds are explicit flags, never silent. Back-projection gives bearing
vectors for future filters.

## 3. PnP: Gauss-Newton on SO(3) × R³

Cost = Σ reprojection². Parameterization: world-frame axis-angle + world
translation; Jacobians analytic (`dq/dφ = −[q]×`, `dq/dt = −R`, pinhole chain
through distortion). Three solver bugs were really one education:

1. **Residual-sign**: `δ = +(J'J)⁻¹J'r` for `r = obs − pred` (negation climbs).
2. **Scale-blind damping**: fixed λ = 1e-6 vs normal entries ~1e6 is no
   damping at all → scaled LM (`λ·diag`) + backtracking line search.
3. **Outlier yank**: one 64 px outlier drags least squares into the wrong
   basin → Huber-IRLS pass 1, hard-gate + unweighted refit, invalid-with-reason
   on starvation (`inlier_starvation`, `rms_above_tolerance`, ...).
4. **Stopping rule**: cost-difference alone stalls under damping → require
   step-norm < 1e-9 too; adaptive λ (down on success for quadratic finish).

## 4. Depth Dilution and Planar Ambiguity (Physics, Not Bugs)

Two findings that looked like bugs: (a) 0.7 m error at 15 m under 0.5 px
noise follows `δz ≈ z²/(f·b_eff)·δpx` with the 1 m plate's weak lever
(`b_eff ≈ 0.2 m`) — subpixel centroiding (0.25 px) halves it; (b) at 20 m the
0.25 m standoff subtends 0.6 px ≈ noise → near-planar depth-flip ambiguity
returns. Rules: standoff-pixels must exceed noise-pixels; acquisition starts
at 15 m. Optical pose NEVER feeds control in M22 (estimates-only discipline
extends to vision: truth comparison is harness-side).

## 5. Verification

11 cases / 90 assertions: principal-point + analytic pixels, flag exclusivity,
distortion direction, bearing round-trip, plate geometry, exact noiseless
recovery (1e-6), noisy bounds, outlier rejection (5/6 inliers, exact refit),
starvation-invalid, monotone sensitivity, config rejection. Oracle (3
functions). Demo: 21-frame approach (worst 0.64 m, all valid), sensitivity to
2 px, occlusion/starvation behavior. Two figures.

## What you should now understand

1. Why four frames, and what breaks if two are confused?
2. Why does the GN sign depend on the residual convention?
3. Why is fixed-λ damping meaningless next to `f = 800` px entries?
4. Why Huber first, hard-gate second (not the reverse)?
5. What ratio governs planar ambiguity, and how do you design out of it?
6. Why must PnP report invalid instead of its best guess?
