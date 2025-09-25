library(dplyr)
library(SimDesign)

## load and combine results
## setwd("paper/simulation/")
dir <- "results/"
files <- list.files(dir)
simresall <- lapply(X = files, FUN = function(file) readRDS(paste0(dir, file)))
simdat <- do.call("rbind", lapply(simresall, FUN = function(x) {
    data.frame(x$condition, x$results)
}))
warnings <- do.call("rbind", lapply(simresall, FUN = function(x) {
    x$warnings
}))
errors <- do.call("rbind", lapply(simresall, FUN = function(x) {
    x$errors
}))

## compute per-condition performance measures
summaries <- simdat |>
    mutate(
        ## sample size difference between treatment group 1 and average
        ## sample size in remaining groups
        N1diff = n1 - (n - n1)/K
    ) |>
    group_by(n, pc, pt1, pt2, K, pH0, capping_eps, burnin, c, method) |>
    summarise(nsim = n(),
              ## mean convergence
              meanconvergence = mean(converged),
              ## expected number of successes
              ENS = mean(nsuccess),
              ENS_sd = sd(nsuccess),
              ENS_mcse = sd(nsuccess)/sqrt(nsim),
              ## expected proportion of successes
              PS = mean(nsuccess/n),
              PS_mcse = sd(nsuccess/n)/sqrt(nsim),
              ## expected number of extreme randomization probabilities
              EExtreme = mean(nextreme),
              EExtreme_mcse = sd(nextreme)/sqrt(nsim),
              ## expected proportion of extreme randomization probabilities
              PExtreme = mean(nextreme/n),
              PExtreme_mcse = sd(nextreme/n)/sqrt(nsim),
              ## expected sample size difference
              EN1diff = mean(N1diff),
              EN1diffL = quantile(N1diff, 0.025),
              EN1diffU = quantile(N1diff, 0.975),
              ## imbalance metric S_0.1
              S01 = mean(N1diff + 0.1*n < 0),
              ## expected allocations to best treatment
              EBest = mean(randbest),
              EBest_mcse = sd(randbest)/sqrt(nsim),
              ## expected proportion of allocations to best treatment
              PBest = mean(randbest/n),
              PBest_mcse = sd(randbest/n)/sqrt(nsim),
              ## expected allocations to treatment 1
              E1 = mean(n1),
              E1_mcse = sd(n1, na.rm = TRUE)/sqrt(nsim),
              ## expected proprtion of allocations to treatment 1
              PE1 = mean(n1/n),
              PE1_mcse = sd(n1/n)/sqrt(nsim),
              ## bias
              biasRD1 = mean(RD1 - pt1 + pc),
              biasRD1_mcse = sd(RD1 - pt1 + pc)/sqrt(nsim),
              ## rejection rate
              rrRD1 = mean(RDlower1 > 0),
              ## coverage
              covRD1 = mean(RDlower1 < pt1 - pc & RDupper1 > pt1 - pc),
    ) |>
    mutate(S01_mcse = sqrt(S01*(1 - S01)/nsim),
           rrRD1_mcse = sqrt(rrRD1*(1 - rrRD1)/nsim),
           covRD1_mcse = sqrt(covRD1*(1 - covRD1)/nsim)) |>
    ungroup()
saveRDS(summaries, "sim-summaries.rds")
