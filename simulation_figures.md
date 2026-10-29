Conditional Q_PC: simulation figures
================

This document regenerates every simulation figure for the conditional
Q_PC / plasticity manuscript from the raw SLiM simulation output. Each
section below mirrors one of the standalone plotting scripts previously
kept in `symmetric/src/` and `negative_rho/src/` (see `SLiM/`); this
file replaces those four scripts as the single source for the figures.
Running it produces this document plus every figure in `figures/` as a
publication-ready `.eps`.

## Symmetric simulations

Trait optima form a symmetric cline (trait1: -1 to 1, trait2: 1 to -1).
Selection is present on both traits (`env2prob = 0.5`) or neither
(`env2prob = 0`, “Neutral”); `rho` (the genetic correlation between the
two traits) is 0.6, 0.75, or 0.9; `sigma_omega` (the width of the
fitness function) is 0.5 or 1.0.

### Q_PC and conditional Q_PC results (Figure 1)

Standard and conditional Q_PC tests on trait1 and trait2, faceted by
selection regime (Neutral vs. Selected) and `sigma_omega`, colored by
`rho`. The x axis uses a two-level nested label (a bracket spanning each
pair): outer label = Q_PC/Conditional Q_PC, inner label = trait1/trait2.

**Figure:** [`figures/sim_results.eps`](figures/sim_results.eps)

### Reaction norm results

Standard and conditional Q_PC tests on the reaction norm (trait1 -
trait2), faceted by selection regime, colored by `rho`, one figure per
`sigma_omega`. A flat (non-nested) x axis, since there is only one trait
(the reaction norm) here.

**Figures:**
[`figures/reaction_norm_sdW0.5.eps`](figures/reaction_norm_sdW0.5.eps)
(sigma_omega = 0.5),
[`figures/reaction_norm_sdW1.0.eps`](figures/reaction_norm_sdW1.0.eps)
(sigma_omega = 1.0)

## Negative-rho simulations

Same symmetric trait-optima cline, but with negative genetic
correlations (`rho` = -0.6, -0.75, -0.9) and selection present on both
traits throughout (`env2prob = 0.5` only – no neutral condition in this
run).

### Q_PC and conditional Q_PC results

Faceted by `sigma_omega` only (there is no selection-regime facet here),
colored by `rho`, with the same nested trait1/trait2 x Q_PC/Conditional
Q_PC x axis as the symmetric Figure 1.

**Figure:**
[`figures/sim_results_negative_rho.eps`](figures/sim_results_negative_rho.eps)

### Reaction norm results

**Figure:**
[`figures/reaction_norm_negative_rho.eps`](figures/reaction_norm_negative_rho.eps)
