#' Bayesian mixture modeling for cylendrical data
#' @importFrom glue glue
#' @param data a data frame, with one column providing data on the
#' angular variable, and one column on the linear variable.
#' @export
fit_cyl_mix <- function(data,
                        ncomp = 2,
                        overfit = "super-impose",
                        seed = 1,
                        sampling_method = "mcmc",
                        init_opt = NULL,
                        angle_col = "angle",
                        line_col = "line",
                        rand_init = (sampling_method == "mcmc"),
                        show_allocation = FALSE,
                        ...)
{

  if (!overfit %in% c("super-impose", "leave-empty")) {
    stop("overfit must be either \'super-impose\' or \'leave-empty\'")
  }

  if (!sampling_method %in% c("mcmc", "vb")) {
    stop("sampling_method must be either \'mcmc\' or \'vb\'")
  }

  lapply(
    c(angle_col, line_col),
    function(xx) {
      if (!xx %in% names(data)) {
        stop(glue::glue("data does not have a column named \'{xx}\'."))
      }
    }
  )

  data <- data[c(line_col, angle_col)]


  if (show_allocation) {
    stan_mix_model <- stanmodels$mixture_w_allocation
  } else {
    stan_mix_model <- stanmodels$mixture
  }

  if (sampling_method == "mcmc") {
    stan_sampling <- rstan::sampling
  } else {
    stan_sampling <- rstan::vb
  }


  alpha0 <- ifelse(overfit == "super-impose", 5.5, 0.01)

  dat.stan <- list(y = data.matrix(data),
                   n_data = nrow(data),
                   n_groups = ncomp,
                   alpha = rep(alpha0, ncomp))


  init <- init1 <- init_kmeans(data, ncomp)
  opt_init_type <- 'kmeans'


  cat("\nOptimizing with kmeans...\n")

  # get MAP

  opt <- opt1 <- tryCatch(
    rstan::optimizing(stan_mix_model,
                      data = dat.stan,
                      seed = seed,
                      init = init1),
    warning = function(w) w,
    error = function(e) e
  )

  opt1_adj <- list()

  if (all(!is(opt, "error"),
          !is(opt, "warning"),
          unlist(opt["return_code"]) == 0)) {
    opt1_adj <- opt1
  }

  if (!is.null(init_opt)) {
    cat("\nOptimizing with provided init..\n")
    opt2 <- tryCatch(rstan::optimizing(stan_mix_model,
                                       data = dat.stan,
                                       seed = seed,
                                       init = init_opt),
                     warning = function(w) w,
                     error = function(e) e)

    opt1_obj <- ifelse(is.null(opt1_adj$value), -Inf, opt1_adj$value)

    if (all(!is(opt2, "error"), !is(opt2, "warning"),
            opt2$return_code == 0, opt2$value > opt1_obj)) {
      cat("\nProvided init yields better MAP than kmeans.. \n")
      opt <- opt2
      init <- init_opt
      opt_init_type <- 'provided'
    } else if (!is.null(opt1_adj$value)) {
      cat("\nProvided init *DOES NOT* yield better MAP than kmeans.. \n")
    }
  }

  if (!is(opt, "error") & !is(opt, "warning")) {
    # MAP etimation is OK!

    map <- sapply(
      c("mu", "sigma",
        "nu", "kappa",
        "pmix"),
      function(xx) {
        opt$par[grepl(xx, names(opt$par))]
      },
      simplify = FALSE,
      USE.NAMES = TRUE
    )

    stan_init_type <- 'map'

    cat("\nOptimization successful. \n")

  } else {
    cat("\nOptimization unsuccessful. \n")
    map <- init
    stan_init_type <- 'kmeans'
  }

  cat("\nStarting sampling..\n")

  dots <- list(...)

  if (is.null(dots$algorithm) & sampling_method == "vb") {
    dots$algorithm <- 'fullrank'
  }

  if (is.null(dots$chains) & sampling_method == "mcmc") {
    dots$chains <- 4
  }

  stan_input_pars <- list(
    object = stan_mix_model,
    data = dat.stan,
    seed = seed)

  if (!rand_init & sampling_method == "mcmc") {
    stan_input_pars$init <- lapply(seq(dots$chains), function(x) map)
  } else if (!rand_init & sampling_method == "vb") {
    stan_input_pars$init <- map
  }


  fit_stan <- do.call(
    stan_sampling,
    c(stan_input_pars, dots)
  )


  out <- list(init = init,
              ncomp = ncomp,
              data = data,
              map = lapply(map, unname),
              stan_map = opt,
              opt_init_type = opt_init_type,
              opt_success = (stan_init_type == "map"),
              stanmodel = fit_stan,
              init_type = stan_init_type)

  class(out) <- "cylmix"
  out
}
