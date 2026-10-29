library(tidyverse)
library(quaint)
library(parallel)

setwd("/mnt/research/josephslab/Husain/cond_Qpc/negative_rho")
args <- commandArgs(trailingOnly = TRUE)
params <- read_table("src/params261021.txt", show_col_types = FALSE) |>
  mutate(index = seed - (min(seed) - 1))

slim_qpc <- function(seed, focal_trait, cor_trait = NULL, test_pcs, var_pcs, directory) {
  phenos <- read.csv(paste(directory, "phenos/", seed, "_pheno.csv", sep = ""))
  genos <- read.table(paste(directory, "genoMat/", seed, "_genoMat.012", sep = ""),
    row.names = 1,
    header = FALSE
  ) / 2
  eigen_k <- eigen(make_k(as.matrix(genos)))
  phenos[, "reaction_norm"] <- phenos[, focal_trait] - phenos[, cor_trait]
  phenos[, "env_mean"] <- rowMeans(phenos[, c("trait1", "trait2")])
  trait1 <- calc_qpc(
    trait_values = phenos[, focal_trait],
    eigen_k = eigen_k,
    test_pcs = test_pcs,
    var_pcs = var_pcs
  )
  trait1_cond <- calc_cond_qpc(
    focal_trait = phenos[, focal_trait],
    cor_trait = phenos[, cor_trait],
    eigen_k = eigen_k,
    test_pcs = test_pcs,
    var_pcs = var_pcs
  )
  trait2 <- calc_qpc(
    trait_values = phenos[, cor_trait],
    eigen_k = eigen_k,
    test_pcs = test_pcs,
    var_pcs = var_pcs
  )
  trait2_cond <- calc_cond_qpc(
    focal_trait = phenos[, cor_trait],
    cor_trait = phenos[, focal_trait],
    eigen_k = eigen_k,
    test_pcs = test_pcs,
    var_pcs = var_pcs
  )
  rxn <- calc_qpc(
    trait_values = phenos[, "reaction_norm"],
    eigen_k = eigen_k,
    test_pcs = test_pcs,
    var_pcs = var_pcs
  )
  rxn_cond <- calc_cond_qpc(
    focal_trait = phenos[, "reaction_norm"],
    cor_trait = phenos[, "env_mean"],
    eigen_k = eigen_k,
    test_pcs = test_pcs,
    var_pcs = var_pcs
  )
  return(list(
    seed = seed,
    trait1 = trait1,
    trait1_cond = trait1_cond,
    trait2 = trait2,
    trait2_cond = trait2_cond,
    rxn = rxn,
    rxn_cond = rxn_cond,
    eigen_k = eigen_k
  ))
}

qpc_list <- mclapply(
  as.numeric(args[2]):as.numeric(args[3]),
  function(x) {
    slim_qpc(
      seed = x,
      focal_trait = "trait1",
      cor_trait = "trait2",
      test_pcs = c(1, 2),
      var_pcs = 200:399,
      directory = args[1]
    )
  },
  mc.cores = detectCores()
)

save.image(file = "qpc_results_negative_rho.RData")
