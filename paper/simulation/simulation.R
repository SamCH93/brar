library(SimDesign)
library(brar)

## uniform priors for exact method
a <- b <- 1
## correlated standard normal prior for normal method
pm <- 0
psd <- 1
rho <- 0.5

## maximum number of treatment groups considered
Kmax <- 3

## fixed objects for SimDesign
fixed_objects <- list("a" = a, "b" = b, "pm" = pm, "psd" = psd, "rho" = rho,
                      "Kmax" = Kmax)

## create fully factorial simulation design
DesignFull <- createDesign(
    ## sample size (low and high power)
    n = c(200, 654),
    ## probability in control group
    pc  = c(0.25),
    ## probability in first treatment group
    pt1 = c(0.25, 0.35, 0.45),
    ## probabilities in remaining treatment groups (if there are any)
    pt2 = c(0.3),
    ## number of treatment groups
    K = seq(1, Kmax),
    ## probability of point null hypothesis
    pH0 = seq(0, 1, 0.25),
    ## threshold 0.5 +- eps where randomization probabilities are capped
    ## (0.5: no capping, 0.4: capping at 0.1 and 0.9)
    capping_eps = c(0.5, 0.4),
    ## length of burn-in (starting phase with equal randomization)
    burnin = c(0, 50),
    ## power-shrinkage (a value of 1 corresponds to no shrinkage)
    c = c("1", "1/2", "i/(2n)"),
    ## exact or normal approximation
    method = c("exact", "normal")
)
## remove combinations that make no sense
Design <- subset(DesignFull,
                 ## "normal" or "exact" produce same result for equal randomization
                 !(pH0 == 1 & method == "normal") &
                 ## power shrinkage only for Thompson sampling (pH0 = 0) settings
                 !(pH0 > 0 & c != "1") &
                 ## capping only for Thompson sampling (pH0 = 0)
                 !(pH0 > 0 & capping_eps != 0.5) &
                 ## burn-in has no effect on equal randomization (pH0 = 1)
                 !(pH0 == 1 & burnin > 0))

## data generation function
Generate <- function(condition, fixed_objects) {
    ## number of groups (1 control + K treatments)
    ngroups <- 1 + condition$K

    ## randomization probability capping limits
    lowercap <- 0.5 - condition$capping_eps
    uppercap <- 0.5 + condition$capping_eps

    ## priors
    if (condition$method == "exact") {
        a1 <- rep(a, ngroups)
        b1 <- rep(b, ngroups)
    }
    if (condition$method == "normal") {
        pm1 <- rep(pm, ngroups - 1)
        psigma1 <- matrix(rho, nrow = ngroups - 1, ncol = ngroups - 1)
        diag(psigma1) <- 1
        psigma1 <- psigma1*psd^2
    }

    ## true outcome probabilities (first is control, second is most effective
    ## treatment, remaining are remaining treatments)
    p <- c(condition$pc, condition$pt1, rep(condition$pt2, condition$K - 1))

    ## create empty vectors to store group allocation, outcomes, and
    ## randomization probability
    group <- integer(condition$n)
    y <- integer(condition$n)
    prand <- matrix(nrow = condition$n, ncol = ngroups)
    colnames(prand) <- paste0("prand", seq(1, ngroups))
    note <- rep(NA_character_, condition$n)

    ## randomize n participants to groups and simulate their outcomes
    for (i in seq_len(condition$n)) {
        if (i <= condition$burnin || condition$pH0 == 1) {
            ## equal randomization if burn-in phase or equal randomization
            prand[i,] <- rep(1/ngroups, ngroups)
        } else {
            ## extract current number patients and events per group
            ycurrent <- sapply(X = seq(1, ngroups), FUN = function(k) {
                sum(y[1:i][which(group[1:i] == k)])
            })
            ncurrent <- sapply(X = seq(1, ngroups), FUN = function(k) {
                length(y[1:i][which(group[1:i] == k)])
            })
            if (condition$method == "exact") {
                res <- try(brar::brar_binomial(y = ycurrent, n = ncurrent,
                                               a0 = a, b0 = b, a = a1, b = b1,
                                               pH0 = condition$pH0),
                           silent = TRUE)
            } else if (condition$method == "normal") {
                dat_glm <- data.frame(y = ycurrent, n = ncurrent,
                                      group = factor(seq(1, ngroups)))
                res <- try({
                    suppressWarnings({
                        fit <- glm(cbind(y, n - y) ~ group, data = dat_glm,
                                   family = binomial)
                    })
                    est <- fit$coef[-1]
                    sigma <- vcov(fit)[-1,-1]
                    brar::brar_normal(estimate = est, sigma = sigma, pm = pm1,
                                      psigma = psigma1, pH0 = condition$pH0)
                }, silent = TRUE)
            } else {
                stop(paste("undefined method:", condition$method))
            }
            ## in case of numerical issues, apply equal randomization
            if (inherits(res, "try-error")) {
                suppressWarnings({
                    res$prand <- rep(1/ngroups, ngroups)
                })
                note[i] <- "applied equal randomization due to numerical issues"
            }
            ## apply power-shrinkage to randomization probabilities
            if (condition$c == "1") {
                c <- 1
            } else if (condition$c == "1/2") {
                c <- 0.5
            } else if (condition$c == "i/(2n)") {
                c <- i/(2*condition$n)
            }
            prand_shrunken <- res$prand^c/sum(res$prand^c)
            ## apply capping to randomization probabilities
            lowcap <- which(prand_shrunken < lowercap)
            upcap <- which(prand_shrunken > uppercap)
            if (length(lowcap) + length(upcap) == 0) {
                prand[i,] <- prand_shrunken
            } else {
                prand_capped <- prand_shrunken
                prand_capped[lowcap] <- lowercap
                prand_capped[upcap] <- uppercap
                ## prand_renorm <- prand_capped/sum(prand_capped)
                ## the above normalization produces probabilities that can be again
                ## below the lower capping range
                ## instead only re-normalize the probabilities above lower limit
                nonlowcap <- which(prand_shrunken > lowercap)
                prand_renorm <- prand_capped
                prand_renorm[nonlowcap] <- prand_capped[nonlowcap]*
                    (1 - sum(prand_capped[lowcap]))/sum(prand_capped[nonlowcap])
                ## UGLY HACK: In rare cases, the re-normalization can still produce
                ## extreme probabilities, so apply a second re-normalization
                ## (using double indexing because of floating point problems
                ## when indexing via prand_renorm < lowercap)
                lowcap2 <- which(prand_renorm[nonlowcap] < lowercap)
                if (length(lowcap2) > 0) {
                    nonlowcap2 <- which(prand_renorm[nonlowcap] > lowercap)
                    prand_renorm[nonlowcap][lowcap2] <- lowercap
                    prand_renorm[nonlowcap][nonlowcap2] <-
                        prand_renorm[nonlowcap][nonlowcap2]*
                        (1 - sum(c(prand_renorm[lowcap],
                                   prand_renorm[nonlowcap][lowcap2])))/
                        sum(prand_renorm[nonlowcap][nonlowcap2])
                }
                prand[i,] <- prand_renorm
            }

        }
        group[i] <- sample(x = seq(1, ngroups), size = 1, replace = TRUE,
                           prob = prand[i,])
        y[i] <- rbinom(n = 1, size = 1, prob = p[group[i]])
    }
    dat <- data.frame(group, y, prand, note)
    return(dat)
}

Analyse <- function(condition, dat, fixed_objects) {
    # proportion of converged cases
    converged <- mean(is.na(dat$note))

    ## number of patients randomized to best group(s)
    truep <- c(condition$pc, condition$pt1, rep(condition$pt2, condition$K - 1))
    bestgroups <- which(truep == max(truep))
    randbest <- sum(dat$group %in% bestgroups)

    ## number of treatment success
    nsuccess <- sum(dat$y)

    ## number of extreme randomization probabilities (either > 0.9 or < 0.1)
    prand <- dat[,paste0("prand", seq(1, condition$K + 1)), drop = FALSE]
    tol <- 1e-10 # numerical tol issues in re-normalization
    nextreme <- sum(sapply(X = seq(1, nrow(prand)), FUN = function(i) {
        any(prand[i,] < 0.1 - tol | prand[i,] > 0.9 + tol)
    }))

    ## risk differences
    K <- condition$K
    n <- sapply(X = seq(1, K + 1),
                FUN = function(i) length(dat$y[dat$group == i]))
    y <- sapply(X = seq(1, K + 1),
                FUN = function(i) sum(dat$y[dat$group == i]))
    p <- y/n
    pSE <- sapply(X = seq(1, K + 1),
                  FUN = function(i) sqrt(p[i]*(1 - p[i])/n[i]))
    RD <- sapply(X = seq(2, K + 1), FUN = function(i) p[i] - p[1])
    RDse <- sapply(X = seq(2, K + 1), FUN = function(i) {
        sqrt(pSE[i]^2 + pSE[1]^2)
    })
    RDp <- 2*(1 - pnorm(abs(RD)/RDse))
    za <- qnorm(p = 0.975)
    RDlower <- RD - za*RDse
    RDupper <- RD + za*RDse

    ## log odds ratios
    a <- y[-1]
    b <- n[-1] - y[-1]
    c <- y[1]
    d <- n[1] - y[1]
    logOR1 <- log(a*d/b/c)
    logORse1 <- sqrt(1/a + 1/b + 1/c + 1/d)
    logOR2 <- log((a + 0.5)*(d + 0.5)/(b + 0.5)/(c + 0.5))
    logORse2 <- sqrt(1/(a + 0.5) + 1/(b + 0.5) + 1/(c + 0.5) + 1/(d + 0.5))
    ## in case of zero cell, use Yates' correction
    logOR <- ifelse(is.finite(logOR1), logOR1, logOR2)
    logORse <- ifelse(is.finite(logORse1), logORse1, logORse2)
    logORp <- 2*(1 - pnorm(abs(logOR)/logORse))
    logORlower <- logOR - za*logORse
    logORupper <- logOR + za*logORse

    ## in case K < Kmax, fill in remaining quantities as NA so that the same
    ## result structure across all simulation conditions
    Kremain <- Kmax - K
    if (Kremain > 0) {
        n <- c(n, rep(NA, Kremain))
        p <- c(p, rep(NA, Kremain))
        RD <- c(RD, rep(NA, Kremain))
        RDse <- c(RDse, rep(NA, Kremain))
        RDp <- c(RDp, rep(NA, Kremain))
        RDlower <- c(RDlower, rep(NA, Kremain))
        RDupper <- c(RDupper, rep(NA, Kremain))
        logOR <- c(logOR, rep(NA, Kremain))
        logORse <- c(logORse, rep(NA, Kremain))
        logORp <- c(logORp, rep(NA, Kremain))
        logORlower <- c(logORlower, rep(NA, Kremain))
        logORupper <- c(logORupper, rep(NA, Kremain))
    }

    ## combine all results
    res <- c(converged,
             randbest, nsuccess, nextreme,
             n, p,
             RD, RDse, RDp, RDlower, RDupper,
             logOR, logORse, logORp, logORlower, logORupper)
    names(res) <- c("converged", "randbest", "nsuccess", "nextreme",
                    paste0("n", seq(0, Kmax)), paste0("p", seq(0, Kmax)),
                    paste0("RD", seq(1, Kmax)), paste0("RDse", seq(1, Kmax)),
                    paste0("RDp", seq(1, Kmax)),
                    paste0("RDlower", seq(1, Kmax)),
                    paste0("RDupper", seq(1, Kmax)),
                    paste0("logOR", seq(1, Kmax)),
                    paste0("logORse", seq(1, Kmax)),
                    paste0("logORp", seq(1, Kmax)),
                    paste0("logORlower", seq(1, Kmax)),
                    paste0("logORupper", seq(1, Kmax)))
    return(res)
}

## ## uncomment to test
## set.seed(42); (condition <- Design[432,]); dat <- Generate(condition); (res <- Analyse(condition, dat))

## save sessionInfo before starting the simulation study
si <- sessionInfo()
saveRDS(si, "sessioninfo-server.rds")

set.seed(143)
nsim <- 10000
ncores <- 100
simres <- runSimulation(
    design = Design,
    generate = Generate,
    analyse = Analyse,
    summarise = NA,
    replications = nsim,
    store_results = TRUE,
    save = TRUE,
    save_results = TRUE,
    filename = "brar-sim.rds",
    control = list(allow_na = TRUE, allow_nan = TRUE),
    fixed_objects = fixed_objects,
    packages = c("brar"),
    parallel = TRUE,
    ncores = ncores
)

## ## load results with
## simres <- readRDS("brar-sim.rds")
## simresall <- SimResults(simres)

## ## get sessionInfo
## attributes(simres)$extra_info$sessionInfo
