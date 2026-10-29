# SLiM

Scripts used to run and analyze the conditional Q_PC simulations, split into
`symmetric/` (symmetric trait-optima cline, `rho` = 0.6/0.75/0.9) and
`negative_rho/` (symmetric trait-optima cline, `rho` = -0.6/-0.75/-0.9).
Both subdirectories contain the same script types, run in the order below;
filenames in `negative_rho/` carry a `_negative_rho` suffix.

1. **`gen_params.r`** — builds the parameter grid (population structure,
   trait optima, `rho`, `sdW`, replicate seeds) and writes it to a
   `params*.txt` file, one row per SLiM run.
2. **`batch_slim.sh`** — SLURM array job; each task reads one row of the
   params file and runs:
   - **`qpc_plasticity.slim`** — the SLiM model: a metapopulation with a
     phenotypic optimum cline across two environments, under selection (or
     neutrality) on two traits.
   - **`tree_seq_pyslim.py`** — loads the resulting tree sequence,
     recapitates and overlays mutations with `pyslim`/`msprime`, and writes
     a VCF (then `vcftools` extracts a genotype matrix from it).
3. **`process_results.sh`** — SLURM wrapper that runs:
   - **`calc_qpc.r`** — for each simulation, computes a kinship matrix from
     the genotype matrix and runs the standard and conditional Q_PC tests on
     each trait and the reaction norm, saving the p-values to a `.RData`
     file.
4. **`summarize_results.r`** — prints summary tables (proportion significant
   by test x selection regime x `sigma_omega` x `rho`) from the tidied Q_PC
   results.
