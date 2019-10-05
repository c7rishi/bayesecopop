#' Bayesian mixture modeling for cylendrical data
#' @importFrom glue glue
#' @param data a data frame, with one column providing data on the
#' @param mixture_type either "finite" or "dirichlet
#' @param overfit either "super-impose" or "leave-empty". Ignored if
#' mixture_type = "dirichlet".
#' angular variable, and one column on the linear variable.
#' @export
fit_cyl_mix <- function(data,
                        mixture_type = "finite",
                        overfit = "super-impose",
                        ncomp = ifelse(mixture_type == "finite" &
                                         overfit == "super-impose",
                                       2, 10),
                        seed = 1,
                        sampling_method = "mcmc",
                        init_opt = NULL,
                        alpha0_dirichlet = 1,
                        alpha0_finite = ifelse(overfit == "super-impose",
                                               5.5, 0.001),
                        angle_col = "angle",
                        line_col = "line",
                        rand_init = (sampling_method == "mcmc"),
                        show_allocation = FALSE,
                        ...)
{

  if (!overfit %in% c("super-impose", "leave-empty")) {
    stop("overfit must be either \'super-impose\' or \'leave-empty\'")
  }

  if (!mixture_type %in% c("finite", "dirichlet")) {
    stop("mixture_type must be either \'finite\' or \'dirichlet\'")
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


  if (show_allocation & mixture_type == "finite") {
    stan_mix_model <- stanmodels$mixture_w_allocation
  } else if (show_allocation & mixture_type == "dirichlet") {
    stan_mix_model <- stanmodels$mixture_dirichlet_w_allocation
  } else if (!show_allocation & mixture_type == "finite") {
    stan_mix_model <- stanmodels$mixture
  } else {
    # !show_allocation & mixture_type == "dirichelt"
    stan_mix_model <- stanmodels$mixture_dirichlet
  }


  if (sampling_method == "mcmc") {
    stan_sampling <- rstan::sampling
  } else {
    stan_sampling <- rstan::vb
  }



  dat.stan <- list(y = data.matrix(data),
                   n_data = nrow(data),
                   n_groups = ncomp,
                   alpha = rep(alpha0_finite, ncomp),
                   alpha0 = alpha0_dirichlet)


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

  if (sampling_method == "vb") {
    if (is.null(dots$algorithm)) {
      dots$algorithm <- 'fullrank'
    }

    if (is.null(dots$output_samples)) {
      dots$output_samples <- 10000
    }
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
