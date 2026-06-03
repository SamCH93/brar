library(SimDesign)
library(brar)

## ---- Design (conditions) ----
Design <- createDesign(rho = seq(0, 0.75, 0.25))

## ---- Fixed constants ----
muc <- 0
mut <- c(0.25, 0, -0.25)
sd <- 1
K <- length(mut)
pm <- rep(0, K)
psd <- 1
pH0 <- 0.5
n <- 200
nincr <- 1
burnin <- 2 * (K + 1)

## ---- Generate ----
Generate <- function(condition, fixed_objects) {
    rho <- condition$rho
    with(fixed_objects, {
        ## Prior mean and covariance
        pm <- rep(0, K)
        psigma <- matrix(rho, nrow = K, ncol = K)
        diag(psigma) <- psd^2

        ## Burn-in data
        datC <- data.frame(
            y = rnorm(ceiling(burnin / (K + 1)), mean = muc, sd = sd),
            group = "Control"
        )
        datT <- do.call("rbind", lapply(seq_len(K), function(k) {
            data.frame(
                y = rnorm(ceiling(burnin / (K + 1)), mean = mut[k], sd = sd),
                group = paste("Treatment", k)
            )
        }))
        dat <- rbind(datC, datT)
        dat$n <- seq_len(nrow(dat))
        dat[, c("pc", paste0("p", 1:K))] <- 1 / (K + 1)

        ## RAR data
        nseq <- seq(tail(dat$n, 1) + 1, n, nincr)
        for (ni in nseq) {
            fit <- lm(y ~ group, data = dat, subset = n < ni)
            estimate <- fit$coef[-1]
            sigma <- vcov(fit)[-1, -1]

            brar <- brar_normal(estimate = estimate, sigma = sigma,
                                pm = pm, psigma = psigma, pH0 = pH0)
            prand <- unname(brar$prand)
            names(prand) <- c("pc", paste0("p", 1:K))

            grpindex <- sample(1:(K + 1), size = nincr, replace = TRUE, prob = prand)
            grp <- c("Control", paste("Treatment", 1:K))[grpindex]
            newdat <- rnorm(nincr, mean = c(muc, mut)[grpindex], sd = sd)

            dati <- data.frame(n = seq(ni, ni + nincr - 1), group = grp, y = newdat)
            dat <- rbind(dat, cbind(dati, data.frame(t(prand))))
        }

        dat
    })
}

## ---- Analyse ----
Analyse <- function(condition, dat, fixed_objects) {
    ## Extract randomization probabilities at interim snapshots
    snapshots <- dat[dat$n %in% c(50, 100, 150, 200), ]
    with(snapshots, {
        ## Mean randomisation probability for treatment 1 at each snapshot
        p1_by_n <- tapply(p1, n, mean)
        setNames(as.numeric(p1_by_n), paste0("p1_n", names(p1_by_n)))
    })
}

## ---- Summarise ----
Summarise <- function(condition, results, fixed_objects) {
    ## Compute mean and quantiles for p1 at each interim snapshot
    do.call("c", lapply(colnames(results), function(col) {
        x <- results[, col]
        q <- quantile(x, probs = c(0.25, 0.5, 0.75))
        setNames(
            c(mean(x), q),
            c(paste0(col, "_mean"), paste0(col, "_q", c(25, 50, 75)))
        )
    }))
}

## ---- Run ----
set.seed(4242)
nsim <- 10000
corsimres <- runSimulation(
    design = Design,
    replications = nsim,
    generate = Generate,
    analyse = Analyse,
    summarise = Summarise,
    fixed_objects = list(
        muc = muc,
        mut = mut,
        sd = sd,
        K = K,
        pm = pm,
        psd = psd,
        pH0 = pH0,
        n = n,
        nincr = nincr,
        burnin = burnin
    ),
    parallel = TRUE,
    ncores = 10,
    store_results = TRUE,
    save = TRUE,
    filename = "simulation/cor-sim.rds",
    packages = c("brar")
)

## ## ---- Plot ----
## corsimres <- readRDS("simulation/cor-sim.rds")
## raw <- SimExtract(corsimres, what = "results")

## plot_dat <- raw |>
##     tidyr::pivot_longer(
##         cols = starts_with("p1_n"),
##         names_to = "n",
##         names_prefix = "p1_n",
##         values_to = "p1"
##     ) |>
##     dplyr::mutate(n = as.integer(n))

## library(ggplot2)
## plot_dat |>
##     dplyr::mutate(
##         nfac = factor(n),
##         rhofac = factor(rho)
##     ) |>
##     ggplot(aes(x = nfac, y = p1, fill = rhofac)) +
##     geom_boxplot(alpha = 0.3) +
##     labs(
##         x = "Total sample size",
##         y = "Randomization probability (group 1)",
##         fill = bquote("Prior correlation" ~ rho)
##     ) +
##     scale_fill_viridis_d() +
##     scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
##     theme_bw() +
##     theme(panel.grid.minor = element_blank(), legend.position = "top")
