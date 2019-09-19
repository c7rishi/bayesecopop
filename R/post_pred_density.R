#' @export
post_pred_density <- function(object, data,
                              angle_col = "angle",
                              line_col = "line",
                              post_sample_size = NULL,
                              angle.max = 2*pi,
                              ...) {
  if (!is(object, "cylmix")) {
    stop("object is not of class \'cylmix\'")
  }

  lapply(
    c(angle_col, line_col),
    function(xx) {
      if (!xx %in% names(data)) {
        stop(glue::glue("data does not have a column named \'{xx}\'."))
      }
    }
  )
  data <- data.matrix(data[c(line_col, angle_col)])
  parnames <- c("mu", "sigma", "nu", "kappa", "pmix")
  samps <- rstan::extract(object$stanmodel)[parnames]
  n_post_samps_full <- length(samps$mu)/object$ncomp


  if (!is.null(post_sample_size)) {
    if (any(post_sample_size <= 0,
            post_sample_size > n_post_samps_full)) {
      stop(glue("\'post_sample_size\' must be a positive number \\
              less than or equal to the total number of \\
              posterior samples in \'object\' ({n_post_samps_full})"))
    }
  }

  if (is.null(post_sample_size)) {
    post_sample_size <- n_post_samps_full
  } else {
    post_sample_size <- ceiling(post_sample_size)
  }


  if (post_sample_size < n_post_samps_full) {
    idx <- sample(n_post_samps_full, post_sample_size)
    samps <- lapply(samps, function(x) x[idx, , drop = FALSE])
  }

  do.call(.d_norm_vm_mix_ave,
          c(list(y = data), samps))
}
