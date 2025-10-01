library(tinytest)
library(brar)

## K = 1 treatment group
estimate <- 0.5
sigma <- 0.1
pm <- 0
psigma <- 1
res <- brar_normal(estimate = estimate, sigma = sigma, pm = pm, psigma = psigma,
                   pH0 = 0.5)

expect_inherits(res, "brar")
expect_true(all(c("data", "prior", "BF_ij", "posterior", "prand") %in% names(res)))
expect_equal(sum(res$posterior), 1, tolerance = 1e-12)
expect_equal(length(res$prand), 2)
expect_equal(names(res$prand), c("Control", "Treatment 1"))
expect_true(all(res$prand >= 0 & res$prand <= 1))

## pH0 = 1 forces equal randomization
res2 <- brar_normal(estimate = 0.3, sigma = 0.2^2, pm = 0, psigma = 1, pH0 = 1)
expect_equal(unname(res2$posterior["H0"]), 1, tolerance = 1e-12)
expect_equal(unname(res2$prand), c(0.5, 0.5), tolerance = 1e-12)

## some invalid inputs
expect_error(
  brar_normal(estimate = c(0.1, 0.2), sigma = matrix(1, 1, 1),
              pm = c(0, 0), psigma = diag(2), pH0 = 0.5)
)
expect_error(
  brar_normal(estimate = NA, sigma = matrix(1, 1, 1),
              pm = 0, psigma = 1, pH0 = 0.5)
)

## K > 1 treatment groups
estimate <- c(0.2, -0.1, 0.3)
sigma <- diag(c(0.05, 0.05, 0.05))
K <- length(estimate)
pm <- rep(0, K)
psigma <- diag(K)
res <- brar_normal(estimate = estimate, sigma = sigma, pm = pm, psigma = psigma,
                   pH0 = 0.5)

## Basic structure
expect_inherits(res, "brar")
expect_true(all(c("data", "prior", "BF_ij", "posterior", "prand") %in% names(res)))

## Posterior probabilities should sum to 1
expect_equal(sum(res$posterior), 1, tolerance = 1e-12)

## Randomization probabilities: Control + 3 treatments
expect_equal(length(res$prand), K + 1)
expect_equal(names(res$prand), c("Control", paste("Treatment", 1:K)))

## All probabilities between 0 and 1
expect_true(all(res$prand >= 0 & res$prand <= 1))

## Bayes factor matrix is square with dimension K + 2
expect_equal(nrow(res$BF_ij), K + 2)
expect_equal(ncol(res$BF_ij), K + 2)

## Symmetry check: BF_ij * BF_ji == 1 (within tolerance)
for (i in 1:(K+2)) {
  for (j in 1:(K+2)) {
    if (i != j) {
      expect_equal(res$BF_ij[i,j] * res$BF_ij[j,i], 1, tolerance = 1e-8)
    }
  }
}
