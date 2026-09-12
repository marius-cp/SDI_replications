## Replication material for the simulation study

- `sim_mean_parameterized.R` and `sim_quantile_parameterized.R` run the reported Monte Carlo experiments. Their outputs use stable, date-free filenames in `data`.
- `sim_mean_plots.R` and `sim_quantile_plots.R` construct the reported rejection-rate figures from the generated result files.
- `plot_true_components.R` constructs the population-component figures and generates `sim_q_compare.rds` when the master file requests a complete rebuild.
- `sim_quantile_unreported.R` runs the additional finite-sample exercise discussed in the paper and writes the unreported diagnostic figure to `plots/unreported_sim.pdf`.
- `bonferroni_dm_comparison.R` contains the supplementary Bonferroni-versus-DM comparison requested during review.

The recommended entry point and expected run times are documented in `../REPRODUCIBILITY.md`.
