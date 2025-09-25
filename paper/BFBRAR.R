## ----"main-setup", include = FALSE--------------------------------------------
## knitr options
library(knitr)
opts_chunk$set(fig.height = 4,
               echo = FALSE,
               warning = FALSE,
               message = FALSE,
               cache = TRUE,
               eval = TRUE,
               fig.align = "center")

## printed digits
options(scipen = 100000)

## should sessionInfo be printed at the end?
Reproducibility <- TRUE

## packages
library(dplyr)
library(SimDesign)
library(brar)
library(mvtnorm)
library(ggplot2)
library(ggh4x)
library(ggpubr)
library(Ternary)

## colors
col1 <- "#D55E00"
col2 <- "#0072B2"
col3 <- "#009E73"


## ----"ternary-plot", fig.height = 3, fig.width = 3----------------------------
## simplex grid
x <- seq(0, 1, length.out = 20)
df <- na.omit(do.call("rbind", lapply(X = x, FUN = function(pp) {
    do.call("rbind", lapply(X = x, FUN = function(pm) {
        if (pp + pm > 1) {
            pm <- NA
            p0 <- NA
        } else {
            p0 <- 1 - pp - pm
        }
        data.frame(pp = pp, pm = pm, p0 = p0, prand = p0/2 + pp)
    }))
})))

par(mar = c(0, 0, 0, 0))
ncols <- 250
colpal <- hcl.colors(n = ncols, palette = "Blue-Red 2", alpha = 0.99, rev = TRUE)
df$col <- colpal[ceiling(df$prand*(ncols - 1) + 1)]
TernaryPlot(alab = bquote("Pr(" * italic(H["-"]) ~ "|" ~ italic(y) * ")" %->% ""),
            blab = bquote("Pr(" * italic(H[0]) ~ "|" ~ italic(y) * ")"  %->% ""),
            clab = bquote("" %<-% "Pr(" * italic(H["+"]) ~ "|" ~ italic(y) * ")"),
            grid.lines = 5, grid.minor.lines = 1,
            axis.labels = list(paste0(seq(0, 100, 20), "%"),
                               paste0(seq(0, 100, 20), "%"),
                               paste0(seq(0, 100, 20), "%")),
            axis.cex = 0.6, lab.cex = 0.8, padding = 0.1)
TernaryPoints(df[, c("pm", "p0", "pp")], pch = 20, col = df$col, cex = 1)
PlotTools::SpectrumLegend("topleft", palette = colpal,
                          legend = paste0(seq(from = 100, to = 0, length.out = 5), "%"),
                          bty = "n", xpd = NA, title = bquote(pi), cex = 0.7, lwd = 15)


## ----"spike-and-slab", fig.height = 3, fig.width = 4--------------------------
par(mar = c(4, 4, 0.5, 0.1))
tseqn <- seq(-3, 0, length.out = 500)
tseqp <- seq(0, 3, length.out = 500)
tseq <- c(tseqn, tseqp)
tau <- 1
mu <- 0
p <- (1 - pnorm((0 - mu)/tau))
pH0 <- 0.5
pHp <- (1 - pH0)*p
pHm <- (1 - pH0)*(1 - p)
densFun <- function(x) {
    pHm*dnorm(x, mu, tau)/(1 - p) + pHp*dnorm(x, mu, tau)/p
}
dens <- densFun(tseq)
plot(tseq, dens, type = "n", xlab = bquote(theta ~ "(treatment effect)"), ylab = "Density",
     ylim = c(0, 1.1*max(c(dens, pH0))), las = 1)
polygon(x = c(tseqn, 0, min(tseqn)),
        y = c(densFun(tseqn), 0, densFun(min(tseqn))), density = 20, angle = 45,
        border = NA, col = adjustcolor(col = col1, alpha = 0.3), lty = 2, lwd = 1)
polygon(x = c(0, tseqp, max(tseqn)), y = c(0, densFun(tseqp), 0), density = 20,
        angle = 135, border = NA, col = adjustcolor(col = col2, alpha = 0.3),
        lty = 2, lwd = 1)
lines(tseqn, densFun(tseqn), col = col1, lwd = 1.5)
lines(tseqp, densFun(tseqp), col = col2, lwd = 1.5)
arrows(x0 = 0, y0 = 0, y1 = pH0, lwd = 1.5, length = 0.06)
## segments(x0 = 0, y0 = 0, y1 = densFun(0), lty = 2)
text(x = 0, y = pH0*1.06, labels = bquote(italic(H)[0]), cex = 1.25)
text(x = 1, y = 0.1, labels = bquote(italic(H)["+"]), cex = 1.25, col = col2)
text(x = -1, y = 0.1, labels = bquote(italic(H)["-"]), cex = 1.25, col = col1)


## ----"check-marginal-likelihood-posterior", eval = FALSE----------------------
# ## check marginal likelihood
# t <- 3
# se <- 2
# mu <- 0
# tau <- 0.1
# taup <- 1/sqrt(1/se^2 + 1/tau^2)
# mup <- (t/se^2 + mu/tau^2)*taup^2
# integrate(f = function(theta) {
#     dnorm(t, theta, se)*dnorm(theta, mu, tau)/pnorm(0, mu, tau, lower.tail = FALSE)
# }, lower = 0, upper = Inf)$value
# dnorm(t, mu, sqrt(se^2 + tau^2))*pnorm(0, mup, taup,lower.tail = FALSE)/
#     pnorm(0, mu, tau, lower.tail = FALSE)
# 
# ## check posteriors
# pH0 <- 0.5
# res <- brar_normal(estimate = t, sigma = se^2, pm = mu, psigma = tau^2,
#                    pH0 = pH0)
# ## closed-form 1 for H+
# 1/(1 + exp(-0.5*(t^2/se^2 - (t - mu)^2/(se^2 + tau^2)))*sqrt(1 + tau^2/se^2)/
#    pnorm(mup/taup)*pH0/(1 - pH0) + pnorm(-mup/taup)/pnorm(mup/taup))
# ## closed-form 2  for H+ (combining quadratic forms)
# 1/(1 + exp(-0.5*((t + mu*se^2/tau^2)^2/(se^2*(1 + se^2/tau^2))) - mu^2/tau^2)*sqrt(1 + tau^2/se^2)/
#    pnorm(mup/taup)*pH0/(1 - pH0) + pnorm(-mup/taup)/pnorm(mup/taup))
# res$posterior[3] # should be the same
# ## closed-form for H-
# 1/(1 + exp(-0.5*(t^2/se^2 - (t - mu)^2/(se^2 + tau^2)))*sqrt(1 + tau^2/se^2)/
#    pnorm(-mup/taup)*pH0/(1 - pH0) + pnorm(mup/taup)/pnorm(-mup/taup))
# res$posterior[1] # should be the same
# ## closed-form for H0
# 1/(1 + (1 - pH0)/pH0*exp(-0.5*((t - mu)^2/(se^2 + tau^2) - t^2/se^2))/sqrt(1 + tau^2/se^2))
# res$posterior[2] # should be the same


## ----"example-normal", fig.height = 3.5, cache = TRUE-------------------------
## simulate adaptive randomization
set.seed(500)
n <- 200
nseq <- seq_len(n)
sd <- 1 # true standard deviation
muc <- 0 # true mean in control group
mut <- 0.25 # true mean in treatment group
datC <- data.frame(y = rnorm(n = n, mean = muc, sd = sd), group = "Control",
                   time = nseq)
datT <- data.frame(y = rnorm(n = n, mean = mut, sd = sd), group = "Treatment",
                   time = nseq)
dat <- rbind(datC, datT)

## priors
pm <- 0
psd <- 1
pH0seq <- seq(0, 1, 0.25)

## perform BRAR for accumulating data
plotDF <- do.call("rbind", lapply(X = seq(5, n), FUN = function(ni) {
    fit <- lm(y ~ group, data = dat, subset = time <= ni)
    estimate <- summary(fit)$coefficients[2,1]
    se <- summary(fit)$coefficients[2,2]
    do.call("rbind", lapply(X = pH0seq, FUN = function(pH0) {
        brar <- brar_normal(estimate = estimate, sigma = se^2, pm = pm,
                            psigma = psd^2, pH0 = pH0)
        res <- data.frame("time" = ni, "pH0" = pH0, "prand" = brar$prand,
                          "group" = names(brar$prand))
        rownames(res) <- NULL
        return(res)
    }))
}))


pH0lvls <- unique(plotDF$pH0)
pH0labs <- as.character(pH0lvls)
pH0labs[pH0lvls == 0] <- "0 (Thompson sampling)"
pH0labs[pH0lvls == 1] <- "1 (equal randomization)"
plotDF$pH0fac <- factor(plotDF$pH0, ordered = TRUE,
                        levels = pH0lvls, labels = pH0labs)
lvls <- c("Control", "Treatment 1")
labs <- c("Control", "Treatment")
## labs <- c(paste("'Control (' ~ mu == ", muc, "~ ')'"),
##           paste("'Treatment (' ~ mu == ", mut, "~ ')'"))
## labs <- c("Control", paste("'Treatment (' * theta ==", mut, "* ')'"))
plotDF$groupLab <- factor(plotDF$group, levels = lvls, labels = labs)
ggplot(data = plotDF, aes(x = time, y = prand, color = pH0fac)) +
    facet_wrap(~ groupLab, labeller = label_parsed) +
    labs(x = "Sample size per group", y = "Randomization probability",
         color = bquote("Pr(" * italic(H)[0] * ")"), linetype = "Method") +
    ## geom_hline(yintercept = 1/(K + 1), lty = 2, alpha = 0.5) +
    geom_step(alpha = 0.8) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "#00000003"),
          legend.position = "top")



## ----"two-dimensional-prior", fig.height = 3.5, fig.width = 3.5---------------
## compute prior density
ngrid <- 300
tseq <- seq(-3, 3, length.out = ngrid)
mu <- c(0, 0)
tau <- 1
rho <- 0.5
sigma <- tau^2*matrix(c(1, rho,
                        rho, 1), ncol = 2, byrow = TRUE)
pH0 <- 0.5
dens <- matrix(t(apply(X = expand.grid(tseq, tseq), 1, FUN = function(x) {
    ## dnorm(x[1], mu[1], tau)*dnorm(x[2], mu[2], tau)
    dmvnorm(x = x, mean = mu, sigma = sigma)
})), nrow = ngrid)*(1 - pH0)
densn <- densp1 <- densp2 <- dens
densn[tseq > 0,] <- NaN
densn[, tseq > 0] <- NaN
densp1[tseq < 0, ] <- NaN
densp1[row(densp1) < col(densp1)] <- NaN
densp2[, tseq < 0] <- NaN
densp2[row(densp2) > col(densp2)] <- NaN

par(mar = c(4, 4, 1, 1))
plot(tseq, tseq, type = "n", xlab = bquote(theta[1] ~ "(effect of treatment 1)"),
     ylab = bquote(theta[2] ~ "(effect of treatment 2)"), las = 1, asp = 1)
contour(x = tseq, y = tseq, z = densn, col = adjustcolor(col1, alpha.f = 0.5),
        drawlabels = FALSE, add = TRUE, lwd = 1.5)
contour(x = tseq, y = tseq, z = densp1, col = adjustcolor(col2, alpha.f = 0.5),
        add = TRUE, lwd = 1.5)
contour(x = tseq, y = tseq, z = densp2, col = adjustcolor(col3, alpha.f = 0.5),
        drawlabels = FALSE, add = TRUE, lwd = 1.5)
segments(x0 = -10, x1 = 0, y0 = 0, lty = 2, col = "#0000004D")
segments(x0 = 0, y0 = -10, y1 = 0, lty = 2, col = "#0000004D")
segments(x0 = 0, x1 = 10, y0 = 0, y1 = 10, lty = 2, col = "#0000004D")
points(0, 0, pch = 4, cex = 1, lwd = 1.5)
text(x = 0, y = 0.3, labels = bquote(italic(H)[0]), cex = 1.25)
text(x = 2, y = -2, labels = bquote(italic(H)["+1"]), cex = 1.25, col = col2)
text(x = -2, y = 2, labels = bquote(italic(H)["+2"]), cex = 1.25, col = col3)
text(x = -2.25, y = -2, labels = bquote(italic(H)["-"]), cex = 1.25, col = col1)


## ----"check-multivariate", eval = FALSE---------------------------------------
# ## simulate and compute probability of H+1
# set.seed(44)
# mu <- c(0.3, -0.4)
# Sigma <- matrix(c(1, 0.8,
#                   0.8, 1),
#                 byrow = TRUE, ncol = 2)
# x <- rmvnorm(n = 1000000, mean = mu, sigma = Sigma)
# 
# ## analytically
# A12 <- matrix(c(-1, 0,
#                 -1, 1), byrow = TRUE, ncol = 2)
# pmvnorm(lower = -Inf, upper = c(0, 0), mean = as.numeric(A12 %*% mu),
#         sigma = A12 %*% Sigma %*% t(A12), keepAttr = FALSE)
# mean(x[,1] > 0 & x[,1] > x[,2])
# ## pmvnorm(lower = -Inf, upper = c(0, 0), mean = mu, sigma = Sigma,
# ##         keepAttr = FALSE)
# ## mean(x[,1] < 0 & x[,2] < 0)
# 
# ## now in 3 dimensions
# mu <- c(0.3, 0, -1)
# Sigma <- matrix(c(1, 0.2, 0,
#                   0.2, 1, 0,
#                   0, 0, 1),
#                 byrow = TRUE, ncol = 3)
# x <- rmvnorm(n = 10000000, mean = mu, sigma = Sigma)
# A13 <- matrix(c(-1, 0, 0,
#                 -1, 1, 0,
#                 -1, 0, 1), byrow = TRUE, ncol = 3)
# pmvnorm(lower = -Inf, upper = c(0, 0, 0), mean = as.numeric(A13 %*% mu),
#         sigma = A13 %*% Sigma %*% t(A13), keepAttr = FALSE)
# mean(x[,1] > 0 & x[,1] > x[,2] & x[,1] > x[,3])
# 
# ## now in 4 dimensions
# mu <- c(0.3, 0, -1, -0.5)
# Sigma <- matrix(c(1, 0.2, 0, 0,
#                   0.2, 1, 0, 0,
#                   0, 0, 1, 0.5,
#                   0, 0, 0.5, 1),
#                 byrow = TRUE, ncol = 4)
# x <- rmvnorm(n = 10000000, mean = mu, sigma = Sigma)
# A14 <- matrix(c(-1, 0, 0, 0,
#                 -1, 1, 0, 0,
#                 -1, 0, 1, 0,
#                 -1, 0, 0, 1), byrow = TRUE, ncol = 4)
# pmvnorm(lower = -Inf, upper = c(0, 0, 0,0 ), mean = as.numeric(A14 %*% mu),
#         sigma = A14 %*% Sigma %*% t(A14), keepAttr = FALSE)
# mean(x[,1] > 0 & x[,1] > x[,2] & x[,1] > x[,3] & x[,1] > x[,4])
# 
# ## check that ratio of determinants is correct
# Tau <- matrix(c(2, 0.6, 0, 0,
#                 0.6, 3, 0.5, 0,
#                 0, 0.5, 4, 0.8,
#                 0, 0, 0.8, 5),
#               byrow = TRUE, ncol = 4)
# det(Sigma + Tau)/det(Sigma)
# det(diag(4) + Tau %*% solve(Sigma))
# det(diag(4) + solve(Sigma) %*% Tau)


## ----"example-multinormal", fig.height = 4, cache = TRUE----------------------
## simulate adaptive randomization
set.seed(42)
n <- 200
nseq <- seq_len(n)
sd <- 1 # true standard deviation
muc <- 0 # true mean in control group
mus <- c(0.25, 0, -0.25) # true means in treatment groups
K <- length(mus) # number of treatments
datC <- data.frame(y = rnorm(n = n, mean = muc, sd = sd), group = "Control",
                   time = nseq)
datT <- do.call("rbind", lapply(seq(1, K), function(k) {
    data.frame(y = rnorm(n = n, mean = mus[k], sd = sd),
               group = paste("Treatment", k), time = nseq)
}))
dat <- rbind(datC, datT)

## priors
pm <- rep(0, K)
rho <- 0.5 # to have equal prior probabilities
psigma <- matrix(rho, nrow = K, ncol = K)
tau <- 1
diag(psigma) <- tau^2
pH0seq <- seq(0, 1, 0.25)

## perform BRAR for accumulating data
plotDF <- do.call("rbind", lapply(X = seq(5, n), FUN = function(ni) {
    fit <- lm(y ~ group, data = dat, subset = time <= ni)
    estimate <- fit$coef[-1]
    sigma <- vcov(fit)[-1,-1]
    do.call("rbind", lapply(X = pH0seq, FUN = function(pH0) {
        brar <- brar_normal(estimate = estimate, sigma = sigma, pm = pm,
                                 psigma = psigma, pH0 = pH0)
        res <- data.frame("time" = ni, "pH0" = pH0, "prand" = brar$prand,
                          "group" = names(brar$prand))
        rownames(res) <- NULL
        return(res)
    }))
}))


pH0lvls <- unique(plotDF$pH0)
pH0labs <- as.character(pH0lvls)
pH0labs[pH0lvls == 0] <- "0 (Thompson sampling)"
pH0labs[pH0lvls == 1] <- "1 (equal randomization)"
plotDF$pH0fac <- factor(plotDF$pH0, ordered = TRUE,
                        levels = pH0lvls, labels = pH0labs)
lvls <- c("Control", paste("Treatment", seq(1, K)))
## labs <- c(paste("'Control (' ~ mu == ", muc, "~ ')'"),
##           paste("'Treatment",  seq(1, K), "(' * mu == ", mus, "* ')'"))
labs <- c("Control",
          paste("'Treatment",  seq(1, K), "(' * theta[", seq(1, K),
                "] == ", mus - muc, "* ')'"))
plotDF$groupLab <- factor(plotDF$group, levels = lvls, labels = labs)
ggplot(data = plotDF, aes(x = time, y = prand, color = pH0fac)) +
    facet_wrap(~ groupLab, labeller = label_parsed) +
    labs(x = "Sample size per group", y = "Randomization probability",
         color = bquote("Pr(" * italic(H)[0] * ")"), linetype = "Method") +
    ## geom_hline(yintercept = 1/(K + 1), lty = 2, alpha = 0.5) +
    geom_step(alpha = 0.8) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "#00000003"),
          legend.position = "top")



## ----"test-derivations", eval = FALSE-----------------------------------------
# set.seed(42)
# a1 <- 2
# b1 <- 5
# a2 <- 1
# b2 <- 3
# a3 <- 5
# b3 <- 7
# nsim <- 1000000
# X <- cbind(rbeta(nsim, a1, b1), rbeta(nsim, a2, b2), rbeta(nsim, a3, b3))
# mean(X[,1] > X[,2] & X[,1] > X[,3])
# 
# intFun <- function(x1) {
#     dbeta(x1, a1, b1) * pbeta(x1, a2, b2) * pbeta(x1, a3, b3)
# }
# integrate(intFun, 0, 1)$value
# 
# ## test whether calculation of Pr(Xi = max(X)) works for arbitrary dimensions
# set.seed(40)
# K <- 10
# a <- sample(seq(1, 10), K, replace = FALSE)
# b <- sample(seq(1, 10), K, replace = FALSE)
# X <- sapply(X = seq_len(K), FUN = function(k) rbeta(nsim, a[k], b[k]))
# mean(apply(X[,1] > X[,-1], 1, prod))
# intFun. <- function(x1) {
#     dbeta(x1, a[1], b[1]) * prod(pbeta(x1, a[-1], b[-1]))
# }
# intFun <- Vectorize(intFun.)
# integrate(intFun, 0, 1)$value
# ## works, yay!


## ----"example-multibinomial", fig.height = 4, cache = TRUE--------------------
## simulate adaptive randomization
set.seed(42)
n <- 200
nseq <- seq_len(n)
sd <- 1 # true standard deviation
pc <- 0.2 # true mean in control group
pt <- c(0.3, 0.2, 0.1) # true means in treatment groups
K <- length(pt) # number of treatments
datC <- data.frame(y = cumsum(rbinom(n = n, size = 1, prob = pc)), n = nseq,
                   group = "Control")
datT <- do.call("rbind", lapply(seq(1, K), function(k) {
    data.frame(y = cumsum(rbinom(n = n, size = 1, prob = pt[k])), n = nseq,
               group = paste("Treatment", k))
}))
dat <- rbind(datC, datT)

## priors
a0 <- 1
b0 <- 1
a <- rep(1, K + 1)
b <- rep(1, K + 1)
pm <- rep(0, K)
rho <- 0.5 # to have equal prior probabilities
psigma <- matrix(rho, nrow = K, ncol = K)
tau <- 1
diag(psigma) <- tau^2
pH0seq <- seq(0, 1, 0.5)

## perform BRAR for accumulating data
plotDF <- do.call("rbind", lapply(X = seq(5, n), FUN = function(ni) {
    fit <- glm(cbind(y, n - y) ~ group, data = subset(dat, n == ni),
               family = binomial)
    est <- fit$coef[-1]
    sigma <- vcov(fit)[-1,-1]
    do.call("rbind", lapply(X = pH0seq, FUN = function(pH0) {
        brar <- brar_binomial(y = subset(dat, n == ni)$y, n = rep(ni, K + 1),
                              pH0 = pH0)
        brarnorm <- brar_normal(estimate = est, sigma = sigma, pm = pm,
                                psigma = psigma, pH0 = pH0)
        res <- rbind(data.frame("time" = ni, "pH0" = pH0, "prand" = brar$prand,
                                "group" = names(brar$prand), method = "exact"),
                     data.frame("time" = ni, "pH0" = pH0,
                                "prand" = brarnorm$prand,
                                "group" = names(brarnorm$prand),
                                method = "normal approximation"))
        rownames(res) <- NULL
        return(res)
    }))
}))


pH0lvls <- unique(plotDF$pH0)
pH0labs <- as.character(pH0lvls)
pH0labs[pH0lvls == 0] <- "0 (Thompson)"
pH0labs[pH0lvls == 1] <- "1 (equal)"
plotDF$pH0fac <- factor(plotDF$pH0, ordered = TRUE,
                        levels = pH0lvls, labels = pH0labs)
lvls <- c("Control", paste("Treatment", seq(1, K)))
labs <- c(paste("'Control (' * theta == ", pc, "* ')'"),
          paste("'Treatment",  seq(1, K), "(' * theta == ", pt, "* ')'"))
plotDF$groupLab <- factor(plotDF$group, levels = lvls, labels = labs)
ggplot(data = plotDF, aes(x = time, y = prand, color = pH0fac,
                          lty = method)) +
    facet_wrap(~ groupLab, labeller = label_parsed) +
    labs(x = "Sample size per group", y = "Randomization probability",
         color = bquote("Pr(" * italic(H)[0] * ")"), linetype = "Method") +
    ## geom_hline(yintercept = 1/(K + 1), lty = 2, alpha = 0.5) +
    geom_step(alpha = 0.8) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "#00000003"),
          legend.position = "top")



## ----"ECMO-analysis", fig.height = 6.3----------------------------------------
## TODO implement randomized play the winner method as comparison
rpw <- function(y1, y0, a = 1, b = 1, g = 1) {
    balls1 <- a + b*y1
    balls0 <- a + g*y0
    p1 <- balls1/(balls1 + balls0)
    c("Control" = 1 - p1, "Treatment 1" = p1)
}

## ECMO trial
y <- c(1, 0, rep(1, 10)) # 1 is survival, 0 is death
treat <- c("ECMO", "control", rep("ECMO", 10))
pH0 <- 0.5
pH0seq <- seq(0, 1, 0.25)
ecmoDF <- do.call("rbind", lapply(X = seq_along(y), FUN = function(i) {
    if (i == 0) {
        y1 <- 0
        n1 <- 0
        y2 <- 0
        n2 <- 0
    } else {
        y1 <- sum(y[1:i][which(treat[1:i] == "control")])
        n1 <- length(y[1:i][which(treat[1:i] == "control")])
        y2 <- sum(y[1:i][which(treat[1:i] == "ECMO")])
        n2 <- length(y[1:i][which(treat[1:i] == "ECMO")])
    }
    a <- y2 + 0.5
    b <- y1 + 0.5
    c <- i - y2 + 0.5
    d <- i - y1 + 0.5
    logOR <- log(a*d/b/c)
    selogOR <- sqrt(1/a + 1/b + 1/c + 1/d)
    condition <- data.frame(y1, n1, y2, n2, ntotal = n1 + n2)
    resbrar <- do.call("rbind", lapply(X = pH0seq, FUN = function(pH0) {
        res <- brar_binomial(y = c(y1, y2), n = c(n1, n2), pH0 = pH0)
        resnor <- brar_normal(estimate = logOR, sigma = selogOR, psigma = 1,
                              pH0 = pH0)
        res <- rbind(data.frame(condition, pH0 = pH0, t(res$posterior),
                                prand = res$prand[2], method = "Exact"),
                     data.frame(condition, pH0 = pH0, t(resnor$posterior),
                                prand = resnor$prand[2], method = "Normal"))
        rownames(res) <- NULL
        return(res)
    }))
    resrpw <- data.frame(condition, pH0 = NA, "H." = NA, "H0" = NA, "H.1" = NA,
                         prand = rpw(y1 = y2, y0 = y1)[2],
                         method = "RPW (original)")
    rbind(resbrar, resrpw)
}))

ecmoDF$pH0lab <- factor(ecmoDF$pH0, ordered = TRUE,
                        levels = pH0seq,
                        labels = c("0 (Thompson)", pH0seq[-c(1,length(pH0seq))],
                                   "1 (equal)"))
plt1 <- ggplot(data = subset(ecmoDF, !is.na(pH0)), #& method != "Normal"),
       aes(x = ntotal, y = prand, linetype = method,
           shape = method)) +
    labs(x = "Total sample size",
         y = "Probability to randomize to ECMO group",
         color = bquote("Pr(" * italic(H)[0] * ")"),
         shape = "", linetype = "") +
    geom_hline(yintercept = 0.5, lty = 2, alpha = 0.5) +
    geom_step(aes(color = pH0lab), alpha = 0.5) +
    geom_point(aes(color = pH0lab)) +
    geom_step(data = subset(ecmoDF, is.na(pH0)), alpha = 0.8) +
    geom_point(data = subset(ecmoDF, is.na(pH0))) +
    scale_x_continuous(breaks = seq(0, 12)) +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())

plt2 <- ggplot(data = ecmoDF,
               aes(x = ntotal, y = H.1, linetype = method, shape = method)) +
    labs(x = "Total sample size",
         y = bquote("Pr(" * italic(H["+"]) ~ "| data)"),
         color = bquote("Pr(" * italic(H)[0] * ")"),
         shape = "", linetype = "") +
    geom_step(aes(color = pH0lab), alpha = 0.5) +
    geom_point(aes(color = pH0lab)) +
    scale_x_continuous(breaks = seq(0, 12)) +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
ggarrange(plt1, plt2, ncol = 1, common.legend = TRUE, legend = "right")


## ----"posterior-probabilities-ECMO"-------------------------------------------
post0 <- subset(ecmoDF, method == "Exact" & ntotal == 12 & pH0 == 0)$H.1
post75 <- subset(ecmoDF, method == "Exact" & ntotal == 12 & pH0 == 0.75)$H.1


## ----"brar-package-demonstration", echo = TRUE--------------------------------
library(brar) # load package

## observed successes and trials in control and 3 treatment groups
y <- c(10, 9, 14, 13)
n <- c(20, 20, 22, 21)
group <- c("control", paste("treatment", seq(1, 3)))

## conduct exact point null Bayesian RAR
brar_binomial(y = y, n = n,
              ## uniform prior for common probability under H0
              a0 = 1, b0 = 1,
              ## uniform priors for all probabilities
              a = c(1, 1, 1, 1), b = c(1, 1, 1, 1),
              ## prior probability of the null hypothesis
              pH0 = 0.5)

## get data into shape for approximate point null Bayesian RAR
fit <- glm(cbind(y, n - y) ~ group, family = binomial) # fit logistic regression
estimate <- fit$coefficients[-1] # remove intercept to get logOR estimates
sigma <- vcov(fit)[-1,-1] # remove intercept to get logOR covariance matrix

## conduct approximate point null Bayesian RAR
pm <- c(0, 0, 0) # set prior mean to zero
## set 0.5 correlated prior covariance so that equal prior probabilities
psigma <- matrix(0.5, nrow = 3, ncol = 3)
diag(psigma) <- 1
brar_normal(estimate = estimate, sigma = sigma, pm = pm, psigma = psigma,
            pH0 = 0.5)


## -----------------------------------------------------------------------------
simres <- readRDS("simulation/brar-sim.rds")
nsim <- unique(simres$REPLICATIONS)


## -----------------------------------------------------------------------------
si <- readRDS("simulation/sessioninfo-server.rds")


## ----"simulation-data"--------------------------------------------------------
## simulation looking at
## - testing: T1E, Power
## - estimation: Bias, variance
## - ethical/patient: Expected number of treatment successes (ENS) / failures (ENF), proportion of patient allocated to best arm

## load results with
simres <- readRDS("simulation/brar-sim.rds")
nsim <- unique(simres$REPLICATIONS)
summaries <- readRDS("simulation/sim-summaries.rds") |>
    mutate(capping = factor(capping_eps, levels = c(0.5, 0.4),
                            labels = c("none", "[0.1,0.9]")),
           method = ifelse(pH0 == 1, "equal", method)) |>
    mutate(trans = case_when(
        capping_eps == 0.5 & c == "1" ~ "none",
        capping_eps == 0.4 & c == "1" ~ "capping",
        capping_eps == 0.5 & c == "1/2" ~ "c=1/2",
        capping_eps == 0.4 & c == "1/2" ~ "capping and c=1/2",
        capping_eps == 0.5 & c == "i/(2n)" ~ "c=i/(2n)",
        capping_eps == 0.4 & c == "i/(2n)" ~ "capping and c=i/(2n)"
                             ),
        trans = factor(trans, levels = c("none", "capping", "c=1/2",
                                         "capping and c=1/2","c=i/(2n)",
                                         "capping and c=i/(2n)"))) |>
    mutate(method = factor(method, levels = c("exact", "normal", "equal"),
                           labels = c("Exact BRAR", "Normal BRAR",
                                      "Equal Randomization"))) |>
    mutate(Klab = paste0("italic(K) == ", K),
           nlab = paste0("italic(n) == ", n),
           burninlab = paste0("'Burn-in' == ", burnin),
           rd1lab = paste0("'RD'[1] == ", pt1 - pc))

## ## warnings and errors
## SimExtract(simres, what = "warnings")
## SimExtract(simres, what = "errors")


## ## compare to results from Robertson et al. (2023)
## summaries |>
##     filter(K == 1, pc == 0.25, pt1 == 0.35, capping_eps == 0.5, burnin == 0,
##            c == "1", pH0 %in% c(0, 1), method == "exact") |>
##     select(n, pH0, method, EN1diff, EN1diffL, EN1diffU, S01, ENS, ENS_sd) |>
##     arrange(-pH0)

## ## ER (n = 200): ENdiff = 0 (-28, 28), S_01 = 0.069, ENS = 60 (6.4)
## ## -> corresponds to pH0 = 1
## ## ER (n = 654): ENdiff = 0 (-50, 50), S_01 = 0.005, ENS = 196 (11.7)
## ## -> corresponds to pH0 = 1
## ## TS (n = 200): ENdiff = 95 (-182, 190), S_01 = 0.137, ENS = 65 (8.5)
## ## -> rougly corresponds to pH0 = 0
## ## TS (n = 654): ENdiff = 461 (-356, 640), S_01 = 0.042, ENS = 220 (17.5)
## ## -> rougly corresponds to pH0 = 0

## TODO create RADAR/spider/star chart that shows multiple performance measures
## at the same time?

## plot parameters
theme_dashboard <- function() {
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "top",
                 panel.grid.minor = element_blank(),
                 panel.grid.major.x = element_blank(),
                 strip.background = element_rect(fill = "#00000003"))
}
dodge_width <- 0.2
errorbar_width <- 0
errorbar_alpha <- 0.8
cols <- c("Exact BRAR" = "#0072B2",
          "Normal BRAR" = "#D55E00",
          "Equal Randomization" = "#000000")
shapes <- c("none" = 19,
            "capping" = 2,
            "c=1/2" = 0,
            "capping and c=1/2" = 5,
            "c=i/(2n)" = 6,
            "capping and c=i/(2n)" = 1)


## ----"convergence"------------------------------------------------------------
maxNC <- summaries |>
    slice_min(order_by = meanconvergence)


## ----"plot-success-rate", dependson = "simulation-data", fig.height = 8, fig.width = 12----
ggplot(data = summaries,
       aes(x = pH0, y = PS, col = method, shape = trans)) +
    facet_nested(nlab + rd1lab ~ burninlab + Klab, labeller = label_parsed) +
    geom_vline(xintercept = seq(-0.125, 1.25, 0.25), alpha = 0.1) +
    # geom_hline(aes(yintercept = pt1), lty = 2, alpha = 0.5) +
    # geom_hline(aes(yintercept = (pt1 + pc)/2), lty = 2, alpha = 0.5) +
    geom_line(alpha = 0.3, position = position_dodge(width = dodge_width)) +
    ## geom_errorbar(aes(ymin = PS - PS_mcse, ymax = PS + PS_mcse),
    ##               width = errorbar_width, alpha = errorbar_alpha,
    ##               position = position_dodge(width = dodge_width)) +
    geom_point(position = position_dodge(width = dodge_width)) +
    scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(labels = scales::percent) +
    scale_color_manual(values = cols) +
    scale_shape_manual(values = shapes) +
    labs(x = bquote("Pr(" * italic(H)[0] * ")"), y = "Mean rate of successes",
         color = "Method", shape = "Modifications") +
    theme_dashboard()


## ----"plot-extreme-probabilities", dependson = "simulation-data", fig.height = 8, fig.width = 12----
ggplot(data = summaries,
       aes(x = pH0, y = PExtreme, color = method, shape = trans)) +
    facet_nested(nlab + rd1lab ~ burninlab + Klab, labeller = label_parsed) +
    geom_vline(xintercept = seq(-0.125, 1.25, 0.25), alpha = 0.1) +
    geom_line(alpha = 0.3, position = position_dodge(width = dodge_width)) +
    ## geom_errorbar(aes(ymin = PExtreme - PExtreme_mcse,
    ##                   ymax = PExtreme + PExtreme_mcse),
    ##               width = errorbar_width, alpha = errorbar_alpha,
    ##               position = position_dodge(width = dodge_width)) +
    geom_point(position = position_dodge(width = dodge_width)) +
    scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(labels = scales::percent) +
    scale_color_manual(values = cols) +
    scale_shape_manual(values = shapes) +
    labs(x = bquote("Pr(" * italic(H)[0] * ")"),
         y = "Mean proportion of extreme randomization probabilities",
         color = "Method", shape = "Modifications") +
    theme_dashboard()


## ----"plot-negative-imbalance", dependson = "simulation-data", fig.height = 8, fig.width = 12----
ggplot(data = filter(summaries, pt1 > pc),
       aes(x = pH0, y = S01, color = method, shape = trans)) +
    facet_nested(nlab + rd1lab ~ burninlab + Klab, labeller = label_parsed) +
    geom_vline(xintercept = seq(-0.125, 1.25, 0.25), alpha = 0.1) +
    geom_line(alpha = 0.3, position = position_dodge(width = dodge_width)) +
    ## geom_errorbar(aes(ymin = S01 - S01_mcse, ymax = S01 + S01_mcse),
    ##               width = errorbar_width, alpha = errorbar_alpha,
    ##               position = position_dodge(width = dodge_width)) +
    geom_point(position = position_dodge(width = dodge_width)) +
    scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(labels = scales::percent) +
    scale_color_manual(values = cols) +
    scale_shape_manual(values = shapes) +
    labs(x = bquote("Pr(" * italic(H)[0] * ")"),
         y = "Proportion of negative sample size differences more than 10%",
         color = "Method", shape = "Modifications") +
    theme_dashboard()


## ----"plot-bias", dependson = "simulation-data", fig.height = 8, fig.width = 12----
ggplot(data = summaries,
       aes(x = pH0, y = biasRD1, color = method, shape = trans)) +
    facet_nested(nlab + rd1lab ~ burninlab + Klab, labeller = label_parsed) +
    geom_vline(xintercept = seq(-0.125, 1.25, 0.25), alpha = 0.1) +
    geom_hline(yintercept = 0, lty = 2, alpha = 0.3) +
    geom_line(alpha = 0.3, position = position_dodge(width = dodge_width)) +
    ## geom_errorbar(aes(ymin = biasRD1 - biasRD1_mcse, ymax = biasRD1 + biasRD1_mcse),
    ##               width = errorbar_width, alpha = errorbar_alpha,
    ##               position = position_dodge(width = dodge_width)) +
    geom_point(position = position_dodge(width = dodge_width)) +
    scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    scale_color_manual(values = cols) +
    scale_shape_manual(values = shapes) +
    labs(x = bquote("Pr(" * italic(H)[0] * ")"),
         y = bquote("Bias (RD"[1] * ")"),
         shape = "Modifications", color = "Method") +
    theme_dashboard()


## ----"plot-coverage", dependson = "simulation-data", fig.height = 8, fig.width = 12----
ggplot(data = summaries,
       aes(x = pH0, y = covRD1, color = method, shape = trans)) +
    facet_nested(nlab + rd1lab ~ burninlab + Klab, labeller = label_parsed) +
    geom_vline(xintercept = seq(-0.125, 1.25, 0.25), alpha = 0.1) +
    geom_hline(yintercept = 0.95, lty = 2, alpha = 0.5) +
    geom_line(alpha = 0.3, position = position_dodge(width = dodge_width)) +
    ## geom_errorbar(aes(ymin = covRD1 - covRD1_mcse, ymax = covRD1 + covRD1_mcse),
    ##               width = errorbar_width, alpha = errorbar_alpha,
    ##               position = position_dodge(width = dodge_width)) +
    geom_point(position = position_dodge(width = dodge_width)) +
    scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(labels = scales::percent) +
    scale_color_manual(values = cols) +
    scale_shape_manual(values = shapes) +
    labs(x = bquote("Pr(" * italic(H)[0] * ")"),
         y = bquote("95% CI Coverage (RD"[1] * ")"),
         shape = "Modifications", color = "Method") +
    theme_dashboard()


## ----"plot-t1e", dependson = "simulation-data", fig.height = 8, fig.width = 12----
ggplot(data = filter(summaries, pt1 == pc),
       aes(x = pH0, y = rrRD1, color = method, shape = trans)) +
    facet_nested(nlab + rd1lab ~ burninlab + Klab, labeller = label_parsed) +
    geom_vline(xintercept = seq(-0.125, 1.25, 0.25), alpha = 0.1) +
    geom_hline(yintercept = 0.025, lty = 2, alpha = 0.5) +
    geom_line(alpha = 0.3, position = position_dodge(width = dodge_width)) +
    ## geom_errorbar(aes(ymin = rrRD1 - rrRD1_mcse, ymax = rrRD1 + rrRD1_mcse),
    ##               width = errorbar_width, alpha = errorbar_alpha,
    ##               position = position_dodge(width = dodge_width)) +
    geom_point(position = position_dodge(width = dodge_width)) +
    scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(labels = scales::percent) + #, limits = c(0, 1)) +
    expand_limits(y = 0) +
    scale_color_manual(values = cols) +
    scale_shape_manual(values = shapes) +
    labs(x = bquote("Pr(" * italic(H)[0] * ")"),
         y = bquote("Type I error rate (RD"[1] == 0 * ")"),
         shape = "Modifications", color = "Method") +
    theme_dashboard()


## ----"plot-power", dependson = "simulation-data", fig.height = 8, fig.width = 12----
ggplot(data = filter(summaries, pt1 > pc),
       aes(x = pH0, y = rrRD1, color = method, shape = trans)) +
    facet_nested(nlab + rd1lab ~ burninlab + Klab, labeller = label_parsed) +
    geom_vline(xintercept = seq(-0.125, 1.25, 0.25), alpha = 0.1) +
    ## geom_hline(yintercept = 0.8, lty = 2, alpha = 0.5) +
    geom_line(alpha = 0.3, position = position_dodge(width = dodge_width)) +
    ## geom_errorbar(aes(ymin = rrRD1 - rrRD1_mcse, ymax = rrRD1 + rrRD1_mcse),
    ##               width = errorbar_width, alpha = errorbar_alpha,
    ##               position = position_dodge(width = dodge_width)) +
    geom_point(position = position_dodge(width = dodge_width)) +
    scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(labels = scales::percent) +
    scale_color_manual(values = cols) +
    scale_shape_manual(values = shapes) +
    labs(x = bquote("Pr(" * italic(H)[0] * ")"),
         y = bquote("Power (RD"[1] == 0 * ")"),
         shape = "Modifications", color = "Method") +
    theme_dashboard()



## ----"sessionInfo2", echo = Reproducibility, results = Reproducibility--------
cat(paste(Sys.time(), Sys.timezone(), "\n"))
sessionInfo()

