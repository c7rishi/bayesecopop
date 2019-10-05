#' @export
estimate_ncomp <- function(object, return_random_compsize = FALSE) {
  nsamp <- nrow(object$data)
  pmix_samps <- rstan::extract(object$stanmodel, pars = "pmix")$pmix
  rand_comp_sizes <- t(
    apply(
      pmix_samps,
      1,
      function(x){
        rmultinom(1, nsamp, x)
      }
    )
  )
  rand_ncomps <- apply(rand_comp_sizes, 1, function(x) sum(x > 0))
  est_ncomp <- as.integer(names(which.max(table(rand_ncomps))))

  if (return_random_compsize) {
    attr(est_ncomp, "rand_compsize") <- rand_comp_sizes
  }

  est_ncomp
}
