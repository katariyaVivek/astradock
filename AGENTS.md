# AstraDock Engineering Instructions

This repository is both an engineering project and an aerospace learning
environment.

Never implement a new mathematical, astrodynamics, GNC, estimation, control,
computer-vision, or machine-learning concept without first explaining the
physical problem and documenting the governing model.

For major concepts, document:

- the physical problem and why the concept is needed,
- coordinate frames and units,
- state variables, inputs, and outputs,
- governing equations and assumptions,
- the numerical method and known failure modes,
- the verification strategy.

Do not hide foundational physics behind third-party aerospace libraries.
During the learning phases, implement foundational algorithms in AstraDock.
External libraries may later be used for independent validation.

Never introduce machine learning where a deterministic classical method should
be understood first. In particular:

- learn an EKF before neural state estimation,
- learn PID and LQR before reinforcement-learning control,
- learn analytical orbital mechanics before learned trajectory models,
- learn geometric pose estimation before learned end-to-end pose estimation.

All simulations must be deterministic when given a random seed. Use SI units
internally unless an exception is explicitly documented. Every major simulation
must keep truth, measurement, estimate, and command state separate.

Before implementing a major subsystem:

1. Inspect the existing code, relevant lessons, and architecture.
2. Write a concise implementation plan.
3. State the physical objective, frames, units, model, assumptions, and
   verification criteria.
4. Implement a small, coherent change with tests.
5. Configure, compile, and run all relevant tests.
6. Inspect and fix project-caused failures and warnings.
7. Summarize what changed and what was learned.

Do not optimize for maximum code output. Optimize for correctness, testability,
physical clarity, maintainability, and educational value. Avoid unrelated
refactors, premature abstraction, hidden global state, ambiguous frames or
units, and undocumented magic numbers.
