# TriadFlux: phase-compatible instantaneous fluid-flux bounds

Research artifact for **Phase-compatible bounds for instantaneous spectral flux in two-dimensional Euler flow**, Shivam Gupta (2026).

[Read the paper](paper/main.pdf).

An exact six-mode Euler initial condition demonstrates that independently maximizing each Fourier triad can overstate attainable instantaneous flux. A balanced example attains only `1/sqrt(2)` of that bound. For larger finite supports, rational relation packings give independently checkable upper bounds.

## Reproduce

Python 3.12 is recommended. From this directory:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python run_experiments.py
python -m unittest discover -s tests -v
python verify_certificates.py
```

The full experiment suite is small and uses CPU only. Figures are saved to `paper/figures/`, raw results to `results/`. Rebuild the manuscript from the repository root with `cd paper && tectonic main.tex`.

The standalone checker needs only Python's standard library:

```sh
python verify_certificates.py
```

## What is proved

- An integer-kernel criterion for attainability of the independent-triad estimate (implemented through Smith normal form).
- A sharp weighted four-triad reduction, with explicitly attaining phases.
- An exact six-mode, real, divergence-free Euler initial-data example.
- A composable upper bound from feasible coefficient allocations on incompatible integer relations.

The equal-weight cycle optimum is classical frustrated-XY mathematics; the paper credits that literature. The fluid construction and certificate workflow are the focused contribution.

## What was measured

- 600 comparisons of grouped, direct-convolution and physical-space flux evaluations.
- 624 complete small-grid weighted-cycle configurations.
- Twelve declared larger-support conditions, including four with no improvement.
- Twelve unit tests, including independent representations, Smith-form torsion cases and tamper rejection.
- Four short-time Navier–Stokes resolution/step-size checks.

The dense-support numerical optimizations are diagnostics, not global certificates. The time integrations are convergence evidence, not validated PDE trajectories. The finite support is an instantaneous initial-data assumption, not an invariant Galerkin subsystem. No statistical-cascade, unresolved-tail, patentability or deployment-performance claim follows from this artifact.

## Layout

- `src/triad_flux.py`: rational coefficients, sharp cycle bound, integer relations, certificate allocation, exact attainability.
- `src/independent.py`: separate direct convolution and physical-space implementation.
- `verify_certificates.py`: standard-library checker with rational upper enclosure of the radical term.
- `tests/`: mathematical, numerical and certificate checks.
- `results/`: raw inputs/results, exact certificates, verification logs and environment versions.
- `paper/`: manuscript source, figures and PDF.
- `docs/`: literature/novelty audit and reproducibility notes.

The verifier validates the bound for the saved cosine polynomial. The analytic derivation and independent equation checks establish its connection to fluid flux.

## License

Code is MIT licensed. Manuscript and original figures: copyright 2026 Shivam Gupta; shared for research reading and citation. Cited works remain under their respective rights.
