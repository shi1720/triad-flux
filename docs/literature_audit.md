# Literature and novelty audit

Search date: 17 October 2026. This is a scoped prior-art assessment, not proof of universal novelty or a patent search.

## Scope selected

Instantaneous spectral flux for a real, finite-support vorticity field on the two-dimensional periodic torus; phase compatibility at fixed mode amplitudes; constructive bounds verifiable from integer relations.

## Closest sources and boundaries

| Source | Established result or method | Boundary of this paper |
|---|---|---|
| [Kraichnan 1967](https://doi.org/10.1063/1.1762301) | Energy/enstrophy conservation, triadic transfer, dual-cascade phenomenology | No claim of a new dual-cascade theory |
| [Murray & Bustamante 2018](https://arxiv.org/abs/1705.08960) | Burgers triad-phase synchronization, transfer efficiency, fixed-amplitude phase dynamics | Phase importance and phase-only models are prior art |
| [Protas, Kang & Bustamante 2024](https://doi.org/10.1103/PhysRevE.109.055104) | Phase alignment in extreme Burgers flows, flux diagnostics and phase randomization | Our optimization is a fixed-amplitude instantaneous 2D Euler bound, with an exact certificate |
| [Benavides & Bustamante 2026](https://arxiv.org/html/2605.03049v1) | Stochastic triad-phase closure for 2D turbulence and cascade-direction predictions | Our deterministic phase compatibility is not a statistical cascade prediction and does not refute weak pairwise correlations |
| [Harper, Bustamante & Nazarenko 2013](https://arxiv.org/abs/1212.3156) | Cluster incidence matrices and nullspace characterization of quadratic wave invariants | We use integer left-null relations of a flux phase map; incidence-matrix algebra is not new |
| [Diep 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC7514657/) | Frustrated XY square has noncollinear optimum with per-bond factor `1/sqrt(2)` | We explicitly credit this optimum; its Euler realization and composable coefficient allocations are the application |
| [Ballouz & Ouellette 2020](https://doi.org/10.1103/PhysRevFluids.5.034603) | Physical-space geometric constraints and transfer efficiency | Our constraint is shared Fourier phase realizability |
| [Kochkov et al. 2021](https://arxiv.org/abs/2102.01010) | Learned numerical fluid solvers | Motivates potential solver/surrogate testing; no learned solver was trained or evaluated here |

Queries included: Euler spectral flux phase synchronization and maximum transfer; Navier–Stokes phase optimization; phase frustration in Fourier triads; phase compatibility constraints; cluster matrices and invariants; sharp triad flux bounds; frustrated XY square optima.

## Claim that survives the audit

A specialized constructive workflow: aggregate the exact instantaneous Euler flux polynomial, detect incompatible signed transfer targets through integer phase relations, realize a sharp six-mode example, and allocate exact penalties to obtain rationally checkable larger-support bounds. The weighted four-cycle proof and elementary deficit estimate are provided self-contained; their ingredients are not presented as newly discovered general mathematics.

No identical Euler construction or this exact certificate workflow was found in the inspected sources. That observation does not establish absence throughout the literature. Acceptance, citation impact, patentability, and market demand have not been established.

## Existing-author manuscript overlap check

Reviewed the first two pages of every top-level PDF in `/Users/shivamgupta/Downloads/research` available at task start. The closest mathematical manuscript is *Frobenius Lens: Exact Experiment Audits for Periodic Cellular Automata*. It concerns finite-field automaton identifiability, not fluid flux. The previous Bernoulli projection and matrix-free lattice-fusion manuscripts concern other mathematical objects. The phase-named volumetry manuscript concerns slice-origin sampling. No prior fluid-mechanics manuscript was identified in that collection. Prior papers were not altered.

## Practical relevance

The demonstrated use is a small, reproducible manufactured test for scientific solvers and spectral surrogates. The exact compatibility certificate can reject a proposed flux target exceeding a necessary physical envelope. No customer interviews, market-size estimates, startup validation, performance savings, or patent opinion are claimed.
