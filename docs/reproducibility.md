# Reproducibility and verification

The mathematical proof establishes global claims. Floating-point tests are diagnostics and implementation checks.

1. `src/triad_flux.py` enumerates the signed-mode terms of the vorticity equation with `Fraction`, aggregating phases and cancelling zeros before constructing certificates.
2. Four-row discovery checks signed pair sums. A general exact rational left-nullspace supplies additional integer witnesses. A full Smith decomposition implements the separate attainability criterion.
3. The emitted allocations are feasible, independent of whether discovery or packing is optimal.
4. `verify_certificates.py` independently checks integer cancellation, sign parity, rational capacities and penalty masses. It encloses sqrt(2) using integer square roots, so its inequality is not dependent on a floating-point optimizer.
5. Separate direct convolution and real-space multiplication verify fluid-flux convention, factors and signs.
6. Short-time runs use explicit RK4 and strict two-thirds dealiasing. They are not interval-validated trajectories.

Raw numerical optimizer starts are reproducible from the declared seed and deterministic run order; optimizer outputs may differ slightly by numerical-library/hardware versions. The exact coefficient and allocation certificates do not depend on those optimizer outputs. Runtime fields in CSV/log files are expected to change on rerun.

The original exploratory search used increasing small integer modes to find a six-mode support with exactly four noncollinear triads and an odd sign product through a radial cutoff. The displayed support is a discovered constructive example. It is not a random sample or evidence that such examples are typical.

All twelve dense conditions are reported, including unchanged bounds. They use L in {2,3,4}, amplitude exponent p in {0,1}, K=L, and energy/enstrophy observables. Optimization is a lower diagnostic, not a certified lower enclosure. No training, external dataset, paid compute, or human-subject evidence is part of the results.
