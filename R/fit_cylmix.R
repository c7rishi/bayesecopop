#' @importFrom glue glue
#' @export
fit_cyl_mix <- function(data, 
                        ncomp = 2, 
                        overfit = "super-impose",
                        show_allocation = TRUE,
                        seed = 1,
                        sampling_method = "mcmc",
                        init_opt = NULL,
                        angle_col = "angle",
                        line_col = "line",
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
        stop(glue::glue("data does not have a column named {xx}."))
      }
    }
  )
  
  
  if (show_allocation) {
    stan_mix_model <- stanmodels$mixture
  } else {
    stan_mix_model <- stanmodel$mixture_w_allocation
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
  
  cat("\nInitializing parameters...\n")
  
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
    
    mu <- opt$par[grepl("mu", names(opt$par))]
    sigma <- opt$par[grepl("sigma", names(opt$par))]
    nu <- opt$par[grepl("nu", names(opt$par))]
    kappa <- opt$par[grepl("kappa", names(opt$par))]
    pmix <- opt$par[grepl("pmix", names(opt$par))]
    
    
    map <- list(mu = mu,
                sigma = sigma,
                nu = nu,
                kappa = kappa,
                pmix = pmix)
    
    vb_init_type <- 'map'
    
    cat("\nOptimization successful. Initializng sampling with MAP..\n")
    
  } else {
    # MAP estimation problematic. Use kmeans
    cat("\nOptimization unsuccessful. Initializng sampling with inits..\n")
    map <- init
    vb_init_type <- 'kmeans'
  }
  
  cat("\nStarting Variational sampling..\n")
  
  dots <- list(...)
  
  if (is.null(dots$algorithm) & method == "vb") {
    dots$algorithm <- 'fullrank'
  } 
  
  fit_stan <- do.call(
    stan_sampling,
    c(list(object = stan_mix_model, 
           data = dat.stan,
           init = map,
           seed = seed),
      dots))
  
  
  
  out <- list(init = init,
              map = lapply(map, unname),
              stan_map = opt,
              opt_init_type = opt_init_type,
              opt_success = (vb_init_type == "map"),
              stanmodel = fit_stan,
              init_type = vb_init_type)
  
  class(out) <- "cylmix"
  out
}
