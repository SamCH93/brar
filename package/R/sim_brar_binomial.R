#' @title Simulate a trial with Bayesian response-adaptive randomization for
#'     binary outcomes
#'
#' @description This function simulates one trial with binary outcomes, a
#'     control group, and one or more treatment groups. Participants are
#'     allocated with Bayesian response-adaptive randomization. Before each
#'     participant, the randomization probabilities are computed with
#'     \code{\link{brar_binomial}} from the data observed so far. Optionally,
#'     the trial can start with a burn-in phase with equal randomization.
#'
#' @param n Total number of participants
#' @param p Vector of true success probabilities. The first element
#'     corresponds to the control group, and the remaining elements correspond
#'     to the treatment groups
#' @param pH0 Prior probability of the null hypothesis. Defaults to
#'     \code{0.5}. Set to \code{0} to obtain Thompson sampling and to \code{1}
#'     to obtain equal randomization
#' @param a0,b0,a,b Beta prior parameters passed to
#'     \code{\link{brar_binomial}}. Default to uniform priors
#' @param burnin Number of participants at the start of the trial that are
#'     allocated with equal randomization. Defaults to \code{0}
#'
#' @return A data frame with one row per participant and the columns
#'     \describe{
#'       \item{\code{id}}{Participant index.}
#'       \item{\code{group}}{Allocated group (\code{0} is control, \code{k} is
#'         treatment \code{k}).}
#'       \item{\code{y}}{Simulated outcome (\code{1} is success).}
#'       \item{\code{prand0}, \code{prand1}, ...}{Randomization probabilities
#'         used to allocate the participant.}
#'       \item{\code{note}}{\code{NA}, or a note if equal randomization was
#'         used because the randomization probabilities could not be
#'         computed.}
#'     }
#'     The simulation settings are stored in the attribute \code{"settings"}.
#'
#' @author Samuel Pawel
#'
#' @seealso \code{\link{oc_brar_binomial}}, \code{\link{brar_binomial}}
#'
#' @examples
#' ## 1 control and 2 treatment groups
#' set.seed(42)
#' dat <- sim_brar_binomial(n = 60, p = c(0.25, 0.45, 0.3), pH0 = 0.5, burnin = 10)
#' head(dat)
#' table(dat$group)
#' plot(dat$id, dat$prand1, type = "s", ylim = c(0, 1),
#'      xlab = "Participant", ylab = "Randomization probability treatment 1")
#'
#' @export
sim_brar_binomial <- function(n, p, pH0 = 0.5, a0 = 1, b0 = 1,
                              a = rep(1, length(p)), b = rep(1, length(p)),
                              burnin = 0) {

    ## input checks
    stopifnot(
        length(n) == 1,
        is.numeric(n),
        is.finite(n),
        n >= 1,
        n == round(n),

        is.numeric(p),
        length(p) >= 2,
        all(is.finite(p)),
        all(p >= 0), all(p <= 1),

        length(pH0) == 1,
        is.numeric(pH0),
        is.finite(pH0),
        0 <= pH0, pH0 <= 1,

        length(a0) == 1, is.numeric(a0), is.finite(a0), a0 > 0,
        length(b0) == 1, is.numeric(b0), is.finite(b0), b0 > 0,
        length(a) == length(p), is.numeric(a), all(is.finite(a)), all(a > 0),
        length(b) == length(p), is.numeric(b), all(is.finite(b)), all(b > 0),

        length(burnin) == 1,
        is.numeric(burnin),
        is.finite(burnin),
        burnin >= 0
    )

    ngroups <- length(p)
    group <- integer(n)
    y <- integer(n)
    prand <- matrix(NA_real_, nrow = n, ncol = ngroups)
    colnames(prand) <- paste0("prand", seq(0, ngroups - 1))
    note <- rep(NA_character_, n)
    ncur <- rep(0, ngroups) # current number of participants per group
    ycur <- rep(0, ngroups) # current number of successes per group

    for (i in seq_len(n)) {
        prandi <- rep(1/ngroups, ngroups)
        if (i > burnin && pH0 < 1) {
            res <- tryCatch(brar_binomial(y = ycur, n = ncur, a0 = a0, b0 = b0,
                                          a = a, b = b, pH0 = pH0),
                            error = function(e) NULL)
            if (is.null(res) || any(!is.finite(res$prand))) {
                note[i] <- "equal randomization due to computational problems"
            } else {
                prandi <- unname(res$prand)
            }
        }
        prand[i,] <- prandi

        ## allocate participant and simulate outcome
        g <- sample.int(n = ngroups, size = 1, prob = prandi)
        yi <- stats::rbinom(n = 1, size = 1, prob = p[g])
        group[i] <- g - 1L
        y[i] <- yi
        ncur[g] <- ncur[g] + 1
        ycur[g] <- ycur[g] + yi
    }

    dat <- data.frame(id = seq_len(n), group = group, y = y, prand,
                      note = note)
    attr(dat, "settings") <- list(outcome = "binomial", n = n, p = p,
                                  pH0 = pH0, a0 = a0, b0 = b0, a = a, b = b,
                                  burnin = burnin)
    return(dat)
}
