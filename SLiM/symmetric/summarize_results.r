library(tidyverse)
setwd("/mnt/research/josephslab/Husain/cond_Qpc/symmetric")

d <- readRDS("figures/sim_results_data.rds")

# --- per test x selection x sigma_omega x rho ---
cat("==== detailed: test x selection x sigma_omega x rho ====\n")
detail <- d |>
  group_by(selection, sigma_omega, test, rho) |>
  summarise(n = n(), n_sig = sum(pval < 0.05), prop_sig = mean(pval < 0.05), .groups = "drop") |>
  arrange(selection, sigma_omega, test, rho)
print(detail, n = 100)

# --- aggregated across rho, per test x selection x sigma_omega ---
cat("\n==== aggregated across rho: test x selection x sigma_omega ====\n")
agg <- d |>
  group_by(selection, sigma_omega, test) |>
  summarise(n = n(), n_sig = sum(pval < 0.05), prop_sig = mean(pval < 0.05), .groups = "drop") |>
  arrange(selection, sigma_omega, test)
print(agg, n = 100)

# --- aggregated across rho AND sigma_omega, per test x selection ---
cat("\n==== aggregated across rho + sigma_omega: test x selection ====\n")
agg2 <- d |>
  group_by(selection, test) |>
  summarise(n = n(), n_sig = sum(pval < 0.05), prop_sig = mean(pval < 0.05), .groups = "drop") |>
  arrange(selection, test)
print(agg2, n = 100)

# --- range of prop_sig across rho within each test x selection x sigma_omega (to describe rho trend) ---
cat("\n==== range across rho ====\n")
rho_range <- detail |>
  group_by(selection, sigma_omega, test) |>
  summarise(min_prop = min(prop_sig), max_prop = max(prop_sig), .groups = "drop")
print(rho_range, n = 100)
