library(tidyverse)
setwd("/mnt/research/josephslab/Husain/cond_Qpc/negative_rho")

d1 <- readRDS("figures/sim_results_negative_rho_data.rds")
d2 <- readRDS("figures/reaction_norm_negative_rho_data.rds")

cat("==== main 4 tests: test x sigma_omega x rho ====\n")
t1 <- d1 |>
  group_by(sigma_omega, test, rho) |>
  summarise(n = n(), n_sig = sum(pval <= 0.05), prop_sig = mean(pval <= 0.05), .groups = "drop") |>
  arrange(sigma_omega, test, rho)
print(t1, n = 100)

cat("\n==== main 4 tests, aggregated across rho: test x sigma_omega ====\n")
t1agg <- d1 |>
  group_by(sigma_omega, test) |>
  summarise(n = n(), n_sig = sum(pval <= 0.05), prop_sig = mean(pval <= 0.05), .groups = "drop")
print(t1agg, n = 100)

cat("\n==== reaction norm tests: test x sigma_omega x rho ====\n")
t2 <- d2 |>
  group_by(sigma_omega, test, rho) |>
  summarise(n = n(), n_sig = sum(pval <= 0.05), prop_sig = mean(pval <= 0.05), .groups = "drop") |>
  arrange(sigma_omega, test, rho)
print(t2, n = 100)

cat("\n==== reaction norm tests, aggregated across rho: test x sigma_omega ====\n")
t2agg <- d2 |>
  group_by(sigma_omega, test) |>
  summarise(n = n(), n_sig = sum(pval <= 0.05), prop_sig = mean(pval <= 0.05), .groups = "drop")
print(t2agg, n = 100)
