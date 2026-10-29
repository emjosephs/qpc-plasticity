library(tidyverse)
setwd("/mnt/research/josephslab/Husain/cond_Qpc/negative_rho/src")

metapopSide <- 10
N <- 50
optimum1 <- "seq(from=-1.0,to=1.0,length=metapopSide)"
optimum2 <- "seq(from=1.0,to=-1.0,length=metapopSide)"
env2prob <- "rep(0.5,metapopSide)" # selected only -- no neutral (env2prob=0) runs
rho <- c(-0.6, -0.75, -0.9)
rep <- 1:100
sdW <- c(0.5, 1.0)

list.params <- list(metapopSide = metapopSide,
                    N = N,
                    optimum1 = optimum1,
                    optimum2 = optimum2,
                    env2prob = env2prob,
                    rho = rho,
                    rep = rep,
                    sdW = sdW)

params <- as_tibble(expand.grid(list.params, stringsAsFactors = F)) |>
  mutate(N_total = (metapopSide^2)*N) |>
  select(-rep)

params$seed <- 2610210000 + (1:nrow(params))

params <- params |>
  select(metapopSide, N, optimum1, optimum2, env2prob, rho, sdW, N_total, seed)

write.table(params,
            file = "params261021.txt",
            quote = F,
            sep = " ",
            row.names = F,
            col.names = T)
