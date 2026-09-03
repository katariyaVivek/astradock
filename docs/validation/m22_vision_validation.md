# M22 Computer Vision & Optical Navigation — Validation Report

## Scope

`cpp/vision/camera_model.hpp` (pinhole, distortion, bearing, fiducials),
`cpp/vision/pnp_pose.hpp` (analytic Jacobians, Huber-LM PnP, gating),
`tests/cpp/test_vision.cpp` (11 cases / 90 assertions), demo
(`tools/vision_demo.cpp`: approach, sensitivity, occlusion), oracle
(`python/audit/independent_vision_reference.py`), telemetry
(`data/m22_*.csv`), figures (`artifacts/figures/m22_*`), lesson
(`docs/lessons/022_vision.md`).

## Analytical checks (all pass)

- Boresight → principal point; (1,−0.5,10) → (400, 200) exact.
- Behind-clip vs out-of-bounds exclusivity; distortion outward + formula.
- Bearing round-trip; world−center−then−rotate; reprojection zero + throw.
- Plate: 1 m edges, 0.25 standoff, <4 rejected.
- Noiseless PnP: exact to 1e-6, range exact, 5/5 inliers.
- Outlier: 5/6 inliers, exact refit. Starvation: invalid + reason.
- Sensitivity monotone; shear-pattern bound documents systematic response.

## Independent verification

Pure-Python oracle (project, bearing, cost); 3 spot checks pass.

## Scenario validation

- Approach 15→5 m @ 0.25 px: 21/21 valid, worst 0.64 m (depth dilution).
- Sensitivity @ 10 m to 2 px: valid throughout, error scales with noise.
- Occlusion: 1-dropped → valid (rms 0.18 px); 2-dropped → correctly refused.

## Development bugs (fixed, kept as findings)

1. GN residual-sign (ascent instead of descent).
2. Scale-blind damping → scaled LM + backtracking.
3. Outlier basin capture → Huber pass 1 + refit.
4. Cost-only stopping → step-norm certificate + adaptive λ.
5. Test-side: wrong behind-camera point; adversarial shear bound; exactness
   expectations under damping.
6. Physics-not-bug × 2: depth dilution lever; planar ambiguity ratio.

## Regression

New `cpp/vision/` + tests/demo only. Full suite 301/301. Ruff clean.

## Limitations

Known correspondences (no association); caller-supplied seed (no global
init); single radial coefficient; ideal point fiducials; pose not yet fused
into navigation (M23+/filter work); no rolling shutter.
