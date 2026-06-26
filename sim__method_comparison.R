setwd("~/Documents/GitHub/qpc-plasticity/")
load('data/Norwich_qpc_inputs.rda')
library(tidyverse)
library(quaint)
library(MBSP)

rho = 0.6
mean_matrix = matrix(0, nrow=140, ncol=2)
var_cov_matrix = matrix(c(1, rho, rho, 1), nrow=2)

###################################################
##                One trait selected             ##
###################################################

run_sim = function(i, rho, mean_matrix, var_cov_matrix, myK, myEig) {
  
  traits = matrix_normal(mean_matrix, myK, var_cov_matrix)
  selection_effect = matrix(c(myEig$vectors[,3]*myEig$values[3],
                              myEig$vectors[,3]*myEig$values[3]*rho), nrow=140)
  selected_traits = traits + selection_effect
  selected_traits = tibble(trait1 = selected_traits[,1], trait2 = selected_traits[,2]) |>
    rowwise() |>
    mutate(rn = trait1 - trait2)
  
  new_traits = tibble(trait1 = selected_traits$trait1) |>
    rowwise() |>
    mutate(trait2 = rho*trait1 + sqrt(1 - rho^2)*rnorm(n = 1, mean = trait1, sd = 1),
           rn = trait1 - trait2)
  
  qpc1     = calc_qpc(trait_values = selected_traits$trait1, eigen_k = myEig, test_pcs = 3, var_pcs = 70:140)
  qpc2     = calc_qpc(trait_values = selected_traits$trait2, eigen_k = myEig, test_pcs = 3, var_pcs = 70:140)
  new_qpc1 = calc_qpc(trait_values = new_traits$trait1,      eigen_k = myEig, test_pcs = 3, var_pcs = 70:140)
  new_qpc2 = calc_qpc(trait_values = new_traits$trait2,      eigen_k = myEig, test_pcs = 3, var_pcs = 70:140)
  qpc_rn     = calc_qpc(trait_values = selected_traits$rn,   eigen_k = myEig, test_pcs = 3, var_pcs = 70:140)
  new_qpc_rn = calc_qpc(trait_values = new_traits$rn,        eigen_k = myEig, test_pcs = 3, var_pcs = 70:140)
  
  cond_qpc1 = calc_cond_qpc(focal_trait = selected_traits$trait1, cor_trait = selected_traits$trait2, 
                            eigen_k = myEig, test_pcs = 3, var_pcs = 70:140)
  cond_qpc2 = calc_cond_qpc(focal_trait = selected_traits$trait2, cor_trait = selected_traits$trait1, 
                            eigen_k = myEig, test_pcs = 3, var_pcs = 70:140)
  
  new_cond_qpc1 = calc_cond_qpc(focal_trait = new_traits$trait1, cor_trait = new_traits$trait2, 
                                eigen_k = myEig, test_pcs = 3, var_pcs = 70:140)
  new_cond_qpc2 = calc_cond_qpc(focal_trait = new_traits$trait2, cor_trait = new_traits$trait1, 
                                eigen_k = myEig, test_pcs = 3, var_pcs = 70:140)
  
  tibble(
    sim = i,
    test = c("qpc1", "qpc2", "new_qpc1", "new_qpc2", "qpc_rn", "new_qpc_rn", 
             "cond_qpc1", "cond_qpc2", "new_cond_qpc1", "new_cond_qpc2"),
    pval = c(qpc1$pvals, qpc2$pvals, new_qpc1$pvals,new_qpc2$pvals, 
             qpc_rn$pvals, new_qpc_rn$pvals, cond_qpc1$pvals, cond_qpc2$pvals, 
             new_cond_qpc1$pvals, new_cond_qpc2$pvals)
  )
}

results = lapply(1:5000, run_sim,
                 rho = rho,
                 mean_matrix = mean_matrix,
                 var_cov_matrix = var_cov_matrix,
                 myK = myK,
                 myEig = myEig) |>
  bind_rows()

results |> 
  mutate(test = factor(test, 
         levels = c("qpc1", "qpc2", "new_qpc1", "new_qpc2", "qpc_rn", "new_qpc_rn", 
                    "cond_qpc1", "cond_qpc2", "new_cond_qpc1", "new_cond_qpc2"))) |>
  group_by(test) |> 
  summarise(reject_rate = mean(pval < 0.05) * 100)

###################################################
##            Reaction norm selected             ##
###################################################

myM = matrix(0, nrow=140, ncol=2)
myV = matrix(c(1, rho, rho, 1), nrow=2)

set.seed(1221)
myZ = matrix_normal(myM, myK, myV) ## sim the neutral traits
myZ = data.frame(myZ)
names(myZ) = c("trait1", "trait2")

#adding selection onto the reaction norm
myZ$rn = myZ[,1] - myZ[,2]
selection_effect = myEig$vectors[,3]*7
myZ$rns = myZ$rn + selection_effect #note need to check direction before adding

calc_qpc(trait_values = myZ$rns, eigen_k = myEig, test_pcs = 3,var_pcs = 70:140)$pvals

# Issue is this is separate from the trait values -- just shows we can add selection onto 
# a reaction norm... can we make it into trait values? Imagine we hold trait 1 steady and 
# add the reaction norm on for trait 2

myZ$trait1 = myZ$trait1 + selection_effect
myZ$rn = myZ$trait1 - myZ$trait2

calc_qpc(trait_values = myZ$rn, eigen_k = myEig, test_pcs = 1:10,var_pcs = 70:140)$pvals[3]

