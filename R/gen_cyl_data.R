#' @export
gen_dat <- function(ndata,
                    parlist = NULL,
                    ncomp = NULL) {
  if (is.null(parlist)) {
    mu <- seq(-2*ncomp, 2*ncomp, length.out = ncomp)
    nu <- BAMBI::minuspi_to_pi(seq(-2*ncomp, 2*ncomp, length.out = 5))
    sigma <- runif(ncomp, 1, 2)
    kappa <- runif(ncomp, 1, 2)
    pmix <- runif(ncomp, 0, 1)
    pmix <- pmix/sum(pmix)
  } else {
    mu <- parlist$mu
    nu <- parlist$nu
    sigma <- parlist$sigma
    kappa <- parlist$kappa
    pmix <- parlist$pmix
    ncomp <- length(pmix)
  }
  gr_indic <- sample(seq_len(ncomp), size = ndata, replace = TRUE,
                     prob = pmix)

  dat <- matrix(NA, nrow = ndata, ncol = 2)

  for (j in seq_len(ncomp)) {
    curr <- which(gr_indic == j)
    ncurr <- length(curr)
    dat[curr, ] <- cbind(rnorm(ncurr, mu[j], sigma[j]),
                         BAMBI::rvm(ncurr, kappa[j], nu[j]))
  }

  colnames(dat) <- c("time", "angle")
  list(dat = as.data.frame(dat),
       pars = list(mu = mu, sigma = sigma,
                   nu = nu, kappa = kappa,
                   pmix = pmix, gr_indic = gr_indic))
}
