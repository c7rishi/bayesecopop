contour_plot_singlepar <- function(data, parlist = NULL, mu = NULL,
                                   sigma = NULL, nu = NULL,
                                   kappa = NULL, pmix = NULL, ...)
{
  colnames_data <- colnames(data)
  if (is.null(colnames_data)) {
    xlab <- ylab <- ""
  } else {
    xlab <- colnames_data[1]
    ylab <- colnames_data[2]
  }

  ncomp <- length(mu)

  if (ncomp > 1) {
    main <- paste("contour plot for fitted", ncomp, "product component mixtures")
  }
  else {
    main <- paste("contour plot for fitted single component model")
  }

  dots <- list(...)

  if (!is.null(dots$main)) {
    main <- dots$main
  }

  npoints <- dots$npoints
  if (is.null(npoints))
    npoints <- 50

  nlevels <- dots$nlevels
  if (is.null(nlevels)) {
    nlevels <- 20
  }

  levels <- dots$levels
  if (is.null(levels)) {
    levels <- exp(seq(-20, 2, length.out = nlevels))
  }

  ypoints <-  seq(min(data[, 2]), max(data[, 2]),
                  length.out = npoints)
  xpoints = seq(min(data[, 1]), max(data[, 1]),
                length.out = npoints)

  coords <- as.matrix(expand.grid(xpoints, ypoints))

  if (!is.null(parlist)) {
    mu <- parlist$mu
    sigma <- parlist$sigma
    nu <- parlist$nu
    kappa <- parlist$kappa
    pmix <- parlist$pmix
  }

  dens <- d_norm_vm_mix(y = coords,
                        mu = mu, sigma = sigma,
                        nu = nu, kappa = kappa,
                        pmix = pmix)
  contour(xpoints, ypoints, matrix(dens, nrow = npoints),
          levels = levels)

  if (is.null(dots$alpha)) {
    dots$alpha <- 0.7
  }

  if (is.null(dots$col)) {
    dots$col <- "red"
  }

  if (is.null(dots$show.data)) {
    dots$show.data <- TRUE
  }



  if (dots$show.data)
    points(data,
           col = scales::alpha(dots$col, dots$alpha),
           pch = 16)

  title(main = main, xlab = xlab, ylab = ylab)
}


d_norm_vm_mix <- function(y, mu, sigma, nu, kappa, pmix) {
  y <- data.matrix(y)
  K <- length(mu)
  dd <- vapply(1:K,
               function(i)
                 log(pmix[i]) +
                 BAMBI::dvm(y[, 2],
                            kappa = kappa[i],
                            mu = nu[i],
                            log = TRUE) +
                 dnorm(y[, 1],
                       mean = mu[i],
                       sd = sigma[i],
                       log = TRUE),
               numeric(nrow(y)))

  max_dd <- max(dd)
  exp(max_dd)*rowSums(exp(dd-max_dd))
}


.d_norm_vm_mix_ave <- function(y, mu, sigma, nu, kappa, pmix) {
  dd <- vapply(1:nrow(mu),
               function(i)
                 d_norm_vm_mix(y, mu[i, ], sigma[i, ],
                               nu[i, ], kappa[i, ], pmix[i, ]),
               numeric(nrow(y)))
  rowMeans(dd)
}


.contour_plot_ave <- function(data, mu, sigma, nu,
                              kappa, pmix, ...)
{
  colnames_data <- colnames(data)
  if (is.null(colnames_data)) {
    xlab <- ylab <- ""
  } else {
    xlab <- colnames_data[1]
    ylab <- colnames_data[2]
  }

  if (is.matrix(mu)) {
    ncomp <- ncol(mu)
  } else {
    ncomp <- 1
  }

  if (ncomp > 1) {
    main <- paste("contour plot for fitted", ncomp, "product component mixtures")
  }
  else {
    main <- paste("contour plot for fitted single component model")
  }

  dots <- list(...)

  if (!is.null(dots$main)) {
    main <- dots$main
  }

  npoints <- dots$npoints
  if (is.null(npoints))
    npoints <- 50

  nlevels <- dots$nlevels
  if (is.null(nlevels)) {
    nlevels <- 20
  }

  levels <- dots$levels
  if (is.null(levels)) {
    levels <- exp(seq(-20, 2, length.out = nlevels))
  }

  ypoints <-  seq(min(data[, 2]), max(data[, 2]),
                  length.out = npoints)
  xpoints = seq(min(data[, 1]), max(data[, 1]),
                length.out = npoints)

  coords <- as.matrix(expand.grid(xpoints, ypoints))
  dens <- .d_norm_vm_mix_ave(y = coords,
                             mu = mu, sigma = sigma,
                             nu = nu, kappa = kappa,
                             pmix = pmix)
  contour(xpoints, ypoints, matrix(dens, nrow = npoints),
          levels = levels)

  if (is.null(dots$alpha)) {
    dots$alpha <- 0.7
  }

  if (is.null(dots$col)) {
    dots$col <- "red"
  }

  if (is.null(dots$show.data)) {
    dots$show.data <- TRUE
  }


  if (dots$show.data)
    points(data,
           col = scales::alpha(dots$col, dots$alpha),
           pch = 16)

  title(main = main, xlab = xlab, ylab = ylab)
}

contour_plot_stan <- function(stanout, data, ...)
{
  samps <- rstan::extract(stanout)
  mu_all <- as.matrix(samps$mu)
  sigma_all <- as.matrix(samps$sigma)
  nu_all <- as.matrix(samps$nu)
  kappa_all <- as.matrix(samps$kappa)
  pmix_all <- as.matrix(samps$pmix)

  .contour_plot_ave(data, mu_all, sigma_all, nu_all,
                    kappa_all, pmix_all, ...)

}


init_kmeans <- function(data, ncomp) {

  data <- data.matrix(data)
  cl <- kmeans(data, centers = ncomp, nstart = 5)
  gr_ids <- lapply(1:ncomp, function(j) which(cl$cluster == j))

  mu <- sigma <- nu <- kappa <- vector("numeric", ncomp)


  for(j in 1:ncomp) {
    if (length(gr_ids[[j]]) == 0) {
      mu[j] <- (j+1) * max(data[, 1]) + runif(1)
      sigma[j] <- 1/100
      nu[j] <- runif(range(data[, 2]))
      kappa[j] <- 100
    } else if (length(gr_ids[[j]]) == 1) {
      mu[j] <- mean(data[gr_ids[[j]], 1])
      sigma[j] <- 1/100
      nu[j] <- circular::mean.circular(data[gr_ids[[j]], 2])
      sigma[j] <- 100
    } else {
      mu[j] <- mean(data[gr_ids[[j]], 1])
      sigma[j] <- sd(data[gr_ids[[j]], 1])

      tmp <- suppressWarnings(circular::mle.vonmises(data[gr_ids[[j]], 2]))
      nu[j] <- as.numeric(tmp$mu)
      kappa[j] <- as.numeric(tmp$kappa)
    }
  }

  n <- vapply(gr_ids, length, 0)
  pmix <- n/sum(n)

  list(mu = mu, sigma = sigma, nu = nu, kappa = kappa, pmix = pmix)

}


est_ncomp_stan <- function(stanout, ...)
{
  samps <- rstan::extract(stanout)
  ncomp <- ncol(samps$mu)
  alloc_all <- samps$allocation
  ncomp_nonempty <- apply(alloc_all, 1,
                          function(x)
                            length(unique(x)))
  as.integer(names(which.max(table(ncomp_nonempty))))
}




contour_plot_gg <- function(parlist,
                            data,
                            ...) {


  `%>%` <- magrittr::`%>%`

  # xvar: linear
  # yvar: angular
  xvar <- names(data)[1]
  yvar <- names(data)[2]

  dots <- list(...)

  angle.max <- dots$angle.max %>% ifelse(is.null(.), 2*pi, .)
  alpha <- dots$alpha %>% ifelse(is.null(.), 0.1, .)
  color <- dots$alpha %>% ifelse(is.null(.), "red", .)
  nvals <- dots$nvals %>% ifelse(is.null(.), 50, .)



  xpoints <- seq(min(data[[xvar]]), max(data[[xvar]]), length.out = nvals)
  ypoints <- seq(min(data[[yvar]]), max(data[[yvar]]), length.out = nvals)
  coords <- as.matrix(expand.grid(xpoints, ypoints * (2*pi)/angle.max))

  # check if MAP or posterior samples

  is_samples <-  "mu[1]" %in% names(parlist)

  if (is_samples) {
    # represent as matrix
    samps <- rstan::extract(parlist)
    parlist_mat <- lapply(samps, as.matrix)
    dens <- .d_norm_vm_mix_ave(y = coords,
                               mu = parlist_mat$mu,
                               sigma = parlist_mat$sigma,
                               nu = parlist_mat$nu,
                               kappa = parlist_mat$kappa,
                               pmix = parlist_mat$pmix)
    ncomp <- dim(parlist_mat$mu)[2]
  } else {
    parlist_mat <- parlist
    dens <- d_norm_vm_mix(y = coords,
                          mu = parlist_mat$mu,
                          sigma = parlist_mat$sigma,
                          nu = parlist_mat$nu,
                          kappa = parlist_mat$kappa,
                          pmix = parlist_mat$pmix)
    ncomp <- length(parlist_mat$mu)
  }

  dat_contour <- tibble::as_tibble(coords) %>%
    dplyr::mutate(dens = dens)


  pl <- ggplot() +
    geom_point(data = data,
               aes_string(x = xvar, y = yvar),
               alpha = alpha, color = color) +
    stat_contour(aes(x = dat_contour$Var1,
                     y = dat_contour$Var2,
                     z = dat_contour$dens,
                     colour = ..level..)) +
    scale_color_gradientn(colours = c("skyblue", "blue", "navy"),
                          name = 'Fitted \n(Post-Predictive) \nDensity') +
    theme_minimal() +
    labs(subtitle = paste(ifelse(is_samples,
                                 "VB",
                                 "MAP"),
                          "fitted",
                          ncomp,
                          "product component mixtures"))
  pl

}


waic.stanfit <- function (x, pars = "log_lik", ..., save_psis = FALSE, cores = getOption("mc.cores",
                                                                                         1))
{
  stopifnot(length(pars) == 1L)
  LLarray <- loo::extract_log_lik(stanfit = x, parameter_name = pars,
                                  merge_chains = FALSE)
  r_eff <- loo::relative_eff(x = exp(LLarray), cores = cores)
  loo::waic.array(LLarray, cores = cores)
}



contour_plot_gg_hexbin <- function(parlist,
                                   data,
                                   ...) {


  `%>%` <- magrittr::`%>%`

  # xvar: linear
  # yvar: angular
  xvar <- names(data)[1]
  yvar <- names(data)[2]


  dots <- list(...)


  angle.max <- dots$angle.max %>% ifelse(is.null(.), 2*pi, .)
  alpha <- dots$alpha %>% ifelse(is.null(.), 0.1, .)
  color <- dots$alpha %>% ifelse(is.null(.), "red", .)
  nvals <- dots$nvals %>% ifelse(is.null(.), 50, .)



  xpoints <- seq(min(data[[xvar]]), max(data[[xvar]]), length.out = nvals)
  ypoints <- seq(min(data[[yvar]]), max(data[[yvar]]), length.out = nvals)
  coords <- as.matrix(expand.grid(xpoints, ypoints))

  # check if MAP or posterior samples

  is_samples <-  "mu[1]" %in% names(parlist)

  if (is_samples) {
    # represent as matrix
    samps <- rstan::extract(parlist)
    parlist_mat <- lapply(samps, as.matrix)
    dens <- .d_norm_vm_mix_ave(y = coords,
                               mu = parlist_mat$mu,
                               sigma = parlist_mat$sigma,
                               nu = parlist_mat$nu,
                               kappa = parlist_mat$kappa,
                               pmix = parlist_mat$pmix)
    ncomp <- dim(parlist_mat$mu)[2]
  } else {
    parlist_mat <- parlist
    dens <- d_norm_vm_mix(y = coords,
                          mu = parlist_mat$mu,
                          sigma = parlist_mat$sigma,
                          nu = parlist_mat$nu,
                          kappa = parlist_mat$kappa,
                          pmix = parlist_mat$pmix)
    ncomp <- length(parlist_mat$mu)
  }

  dat_contour <- tibble::as_tibble(coords) %>%
    dplyr::mutate(Var2 = Var2,
                  dens = dens) %>%
    mutate(Var2 = Var2 * angle.max/(2*pi)) # rescale

  data[[yvar]] <- data[[yvar]] * angle.max/(2*pi) # rescale


  pl <-
    ggplot() +
    geom_hex(data = data,
             aes_string(x = xvar, y = yvar)) +
    stat_contour(aes(x = dat_contour$Var1,
                     y = dat_contour$Var2,
                     z = dat_contour$dens,
                     colour = ..level..)) +
    scale_color_gradientn(colours = c("cyan", "blue", "black"),
                          name = 'Fitted \n(Post-Pred) \nDensity') +
    scale_fill_gradientn(colours = c("white", "red", "chocolate4"),
                         name = 'Data \nCount') +
    theme_minimal() +
    labs(subtitle = paste(ifelse(is_samples,
                                 "VB",
                                 "MAP"),
                          "fitted",
                          ncomp,
                          "product component mixtures")) +
    coord_flip()


  pl

}


surface_plot_plotly <- function(parlist,
                                data,
                                ...) {


  `%>%` <- magrittr::`%>%`

  # xvar: linear
  # yvar: angular
  xvar <- names(data)[1]
  yvar <- names(data)[2]


  dots <- list(...)


  angle.max <- dots$angle.max %>% ifelse(is.null(.), 2*pi, .)
  alpha <- dots$alpha %>% ifelse(is.null(.), 0.1, .)
  color <- dots$alpha %>% ifelse(is.null(.), "red", .)
  nvals <- dots$nvals %>% ifelse(is.null(.), 50, .)



  xpoints <- seq(min(data[[xvar]]), max(data[[xvar]]), length.out = nvals)
  ypoints <- seq(min(data[[yvar]]), max(data[[yvar]]), length.out = nvals)
  coords <- as.matrix(expand.grid(xpoints, ypoints * (2*pi)/angle.max))

  # check if MAP or posterior samples

  is_samples <-  "mu[1]" %in% names(parlist)

  if (is_samples) {
    # represent as matrix
    samps <- rstan::extract(parlist)
    parlist_mat <- lapply(samps, as.matrix)
    dens <- .d_norm_vm_mix_ave(y = coords,
                               mu = parlist_mat$mu,
                               sigma = parlist_mat$sigma,
                               nu = parlist_mat$nu,
                               kappa = parlist_mat$kappa,
                               pmix = parlist_mat$pmix)
    ncomp <- dim(parlist_mat$mu)[2]
  } else {
    parlist_mat <- parlist
    dens <- d_norm_vm_mix(y = coords,
                          mu = parlist_mat$mu,
                          sigma = parlist_mat$sigma,
                          nu = parlist_mat$nu,
                          kappa = parlist_mat$kappa,
                          pmix = parlist_mat$pmix)
    ncomp <- length(parlist_mat$mu)
  }

  dat_contour <- tibble::as_tibble(coords) %>%
    dplyr::mutate(dens = dens)


  dens_mat <- matrix(dens,
                     nrow = length(xpoints),
                     ncol = length(ypoints))


  subtitle <- paste(ifelse(is_samples,
                           "VB",
                           "MAP"),
                    "fitted",
                    ncomp,
                    "product component mixtures")

  title <-  ifelse(is.null(dots$title),
                   "",
                   paste0(dots$title, "\n")) %>%
    paste0(subtitle)


  Density <- dens_mat

  pl <-
    plotly::plot_ly(x = ~xpoints,
                    y = ~ypoints,
                    z = ~Density) %>%
    plotly::add_surface() %>%
    plotly::layout(
      title = title,
      scene = list(xaxis = list(title = xvar),
                   yaxis = list(title = yvar),
                   zaxis = list(title = "Density"))
    )

  pl

}



surface_plot_lattice <- function(parlist,
                                 data,
                                 ...) {


  `%>%` <- magrittr::`%>%`

  # xvar: linear
  # yvar: angular
  xvar <- names(data)[1]
  yvar <- names(data)[2]


  dots <- list(...)


  angle.max <- dots$angle.max %>% ifelse(is.null(.), 2*pi, .)
  alpha <- dots$alpha %>% ifelse(is.null(.), 0.1, .)
  color <- dots$alpha %>% ifelse(is.null(.), "red", .)
  nvals <- dots$nvals %>% ifelse(is.null(.), 50, .)
  xlab <- dots$xlab %>% ifelse(is.null(.), yvar, .)
  ylab <- dots$ylab %>% ifelse(is.null(.), xvar, .)



  xpoints <- seq(min(data[[xvar]]), max(data[[xvar]]), length.out = nvals)
  ypoints <- seq(min(data[[yvar]]), max(data[[yvar]]), length.out = nvals)
  coords <- as.matrix(expand.grid(xpoints, ypoints))

  # check if MAP or posterior samples

  is_samples <-  "mu[1]" %in% names(parlist)

  if (is_samples) {
    # represent as matrix
    samps <- rstan::extract(parlist)
    parlist_mat <- lapply(samps, as.matrix)
    dens <- .d_norm_vm_mix_ave(y = coords,
                               mu = parlist_mat$mu,
                               sigma = parlist_mat$sigma,
                               nu = parlist_mat$nu,
                               kappa = parlist_mat$kappa,
                               pmix = parlist_mat$pmix)
    ncomp <- dim(parlist_mat$mu)[2]
  } else {
    parlist_mat <- parlist
    dens <- d_norm_vm_mix(y = coords,
                          mu = parlist_mat$mu,
                          sigma = parlist_mat$sigma,
                          nu = parlist_mat$nu,
                          kappa = parlist_mat$kappa,
                          pmix = parlist_mat$pmix)
    ncomp <- length(parlist_mat$mu)
  }

  subtitle <- paste(ifelse(is_samples,
                           "VB",
                           "MAP"),
                    "fitted",
                    ncomp,
                    "product component mixtures")

  title <-  ifelse(is.null(dots$title),
                   "",
                   paste0(dots$title, "\n")) %>%
    paste0(subtitle)


  if (is.null(dots$screen)) {
    dots$screen <- list(z = 45, x = -60)
  }

  dots$xlab <- NULL
  dots$ylab <- NULL

  pl <- do.call(BAMBI:::basic_surfaceplot,
                c(list(xpoints = ypoints * angle.max/(2*pi),
                       ypoints = xpoints,
                       denmat = t(dens_mat),
                       xlab = xlab,
                       ylab = ylab,
                       main = title),
                  dots))


  pl

}

# only need data if xpoints and ypoints are missing
# coords = xpoints x ypoints (cartesian product)
# if coords is non NULL, xpoints & ypoints & data are ignored
.eval_density_parlist <- function(parlist,
                                  data = NULL,
                                  xpoints = NULL,
                                  ypoints = NULL,
                                  coords = NULL,
                                  nvals = 50,
                                  angle.max = 2*pi,
                                  ...) {
  dots <- list(...)
  if (is.null(coords)) {
    if (is.null(xpoints)) {
      xvar <- names(data)[1]
      xpoints <- seq(min(data[[xvar]]), max(data[[xvar]]), length.out = nvals)
    }

    if (is.null(ypoints)) {
      yvar <- names(data)[2]
      ypoints <- seq(min(data[[yvar]]), max(data[[yvar]]), length.out = nvals)
    }


    coords <- as.matrix(expand.grid(xpoints, ypoints * (2*pi)/angle.max))
  }
  # check if MAP or posterior samples

  is_samples <-  any(
    "mu[1]" %in% names(parlist),
    is.matrix(parlist$mu),
    length(parlist$mu) > 1
  )

  if (is_samples) {
    # represent as matrix
    if (is(parlist, "stanfit")) {
      samps <- rstan::extract(parlist)
    } else {
      samps <- parlist
    }
    parlist_mat <- lapply(samps, as.matrix)
    dens <- .d_norm_vm_mix_ave(y = coords,
                               mu = parlist_mat$mu,
                               sigma = parlist_mat$sigma,
                               nu = parlist_mat$nu,
                               kappa = parlist_mat$kappa,
                               pmix = parlist_mat$pmix)
    ncomp <- dim(parlist_mat$mu)[2]
  } else {
    parlist_mat <- parlist
    dens <- d_norm_vm_mix(y = coords,
                          mu = parlist_mat$mu,
                          sigma = parlist_mat$sigma,
                          nu = parlist_mat$nu,
                          kappa = parlist_mat$kappa,
                          pmix = parlist_mat$pmix)
    ncomp <- length(parlist_mat$mu)
  }

  if (!is.null(xpoints) & !is.null(ypoints)) {
    dens_mat <- matrix(dens,
                       nrow = nrow(xpoints),
                       ncol = ncol(ypionts))

  } else {
    dens_mat <- dens
  }

  list(denmat = dens_mat, ncomp = ncomp,
       xpoints = xpoints, ypoints = ypoints)
}

overlay_surface_plotly <- function(parlists, data, ...) {
  `%>%` <- magrittr::`%>%`

  dots <- list(...)

  angle.max <- dots$angle.max %>% ifelse(is.null(.), 2*pi, .)
  nvals <- dots$nvals %>% ifelse(is.null(.), 50, .)
  xvar <- names(data)[1]
  yvar <- names(data)[2]
  xpoints <- seq(min(data[[xvar]]), max(data[[xvar]]), length.out = nvals)
  ypoints <- seq(min(data[[yvar]]), max(data[[yvar]]), length.out = nvals)

  # buttons will be created in 3D surfaces
  # based on these names
  allsites <- names(parlists)
  if (any(is.null(allsites))) {
    allsites <- paste("Group", 1:length(parlists))
  }
  nsites <- length(allsites)


  alldens <- lapply(parlists, .eval_density_parlist,
                    xpoints = xpoints, ypoints = ypoints,
                    angle.max = angle.max, nvals = nvals)
  denmats <- lapply(alldens, "[[", "denmat")
  names(alldens) <- names(denmats) <- allsites


  title <-  ifelse(is.null(dots$title),
                   "",
                   dots$title)

  showscale  <-  ifelse(is.null(dots$showscale),
                        FALSE, .)

  # need to eval(parse(...)) to evaluate
  # plotly objects in loop
  pstring <-
    "pl <- plotly::plot_ly(x = ~xpoints,
                           y = ~ypoints,
                           z = ~tmp,
                           showscale = showscale)"

  for (j in 1:nsites) {
    tmp <- denmats[[j]]
    assign(paste0("denmat_", j, sep=""), tmp)

    pstring <- pstring %>%
      glue::glue(" %>% plotly::add_surface(
                    z = ~{eval(paste0(\"denmat_\", j))},
                    name = \"{allsites[j]}\",
                    showscale = showscale
                    )")
  }

  eval(parse(text = pstring))

  button_method <- 'update'

  buttons <- purrr::map2(
    allsites,
    seq(nsites),
    function(x, ii) {
      # indic <- c(TRUE, FALSE)
      indic <- rep(FALSE, nsites)
      indic[ii] <- TRUE
      purrr::map2(c(x, paste0("-", x)),
                  list(indic, !indic),
                  function(y, ind)
                    list(method = button_method,
                         args = list(list(visible = ind)),
                         label = y))
    }
  )

  buttons <- unlist(buttons, recursive = FALSE)
  nbuttons <- length(buttons)

  buttons[[nbuttons+1]] <- list(
    method = button_method,
    args = list(list(visible = rep(TRUE, nsites))),
    label = "All"
  )

  buttons[[nbuttons+2]] <- list(
    method = button_method,
    args = list(list(visible = rep(FALSE, nsites))),
    label = "None"
  )

  xlab <- dots$xlab %>% ifelse(is.null(.), xvar, .)
  ylab <- dots$ylab %>% ifelse(is.null(.), yvar, .)


  pl <- pl %>%
    plotly::layout(
      title = title,
      scene = list(xaxis = list(title = xvar),
                   yaxis = list(title = yvar),
                   zaxis = list(title = "Fitted Density")),
      updatemenus = list(
        list(
          type = "buttons",
          y = 0.8,
          buttons = buttons))
    )

  pl
}


# post predicitive sample from one parameter theta_i
.post_pred_sample_onepar <- function(mu, sigma, nu, kappa, pmix) {
  ncomp <- length(pmix)
  idx <- sample(seq_len(ncomp), 1, prob = pmix)
  line <- rnorm(n = 1, mean = mu[idx], sd = sigma[idx])
  angle <- BAMBI::rvm(n = 1, mu = nu[idx], kappa = kappa[idx])
  list(line = line, angle = angle)
}



# posterior predicitve sampling given MCMC iterations
# of mixture model parameters
# parlist = either stanfit (MCMC/VB simulations) OR
# rstan::extract(stanfit) <-- faster

.post_pred_sample_parlist <- function(n = 1,
                                      parlist) {

  if (is(parlist, "stanfit")) {
    allpars <- rstan::extract(parlist)
  } else {
    allpars <- parlist
  }

  nsim <- nrow(allpars$pmix)


  sample_idx <- sample(seq_len(nsim), n, replace = n > nsim)

  pars_by_idx <- lapply(
    sample_idx,
    function(ii) {
      lapply(allpars[c("mu", "sigma",
                       "nu", "kappa",
                       "pmix")],
             function(x) x[ii, ])
    }
  )
  # produces a list of size n, each element is a further list of
  # one iteration of the parameters


  pp_samples_list <- lapply(
    pars_by_idx,
    function(xx) {
      do.call(.post_pred_sample_onepar, xx)
    }
  )

  r_line <- vapply(pp_samples_list, "[[", 0,  "line")
  r_angle <- vapply(pp_samples_list, "[[", 0,  "angle")

  list(line = r_line, angle = r_angle)
}



.list_by_row <- function(mat) {
  lapply(1:nrow(mat), function(ii) mat[ii, ])
}


# 1. generates pp samples from parlist1 and parlist2
# 2. evaluates density of the pp samples at each param
# in parlist1 and also parlist2
# makes a table with columns: -->
# line, angle, mod, den_indiv_mod_1_1:length(parlist1),
# den_indiv_mod_2_1:length(parlist2)


.table_pp_sample_density <- function(n = 1000,
                                     parlist1,
                                     parlist2) {
  allparlists <- list(`mod1` = parlist1,
                      `mod2` = parlist2) %>%
    lapply(
      function(parlist) {
        if (is(parlist, "stanfit")) {
          samps <- rstan::extract(parlist)
        } else {
          samps <- parlist
        }
        samps
      }
    )


  pp_samps <- purrr::imap(
    allparlists,
    function(parlist, ii) {
      tmp <- do.call(
        cbind,
        .post_pred_sample_parlist(n, parlist)
      )
      tmp_dt <- as.data.frame(tmp)
      tmp_dt$mod <- ii
      tmp_dt
    }
  ) %>%
    do.call(rbind, .)


  den_allpar_mat <- purrr:::imap(
    allparlists,
    function(parlist, ii) {
      tmp <- purrr::pmap(
        parlist[c('mu', 'sigma', 'nu',
                  'kappa', 'pmix')] %>%
          lapply(.list_by_row),
        function(mu, sigma, nu, kappa, pmix) {
          d_norm_vm_mix(
            y = data.matrix(pp_samps[, 1:2]),
            mu = mu,
            sigma = sigma,
            nu = nu,
            kappa = kappa,
            pmix = pmix
          )
        }
      ) %>%
        do.call(cbind, .)
      colnames(tmp) <- paste("den_indiv", ii,
                             seq_len(ncol(tmp)), sep = "_")
      as.data.frame(tmp)
    }
  ) %>%
    do.call(cbind,. )

  cbind(pp_samps, den_allpar_mat)

}




# # # randomly interchange model labels (1 & 2's)
# # .rand_swap_mod_samps <- function(pp_df) {
# #   `%>%` <- magrittr::`%>%`
# #
# #   pp_df <- data.matrix(pp_df)
# #   df_names <- colnames(pp_df)
# #
# #   mod_1_names <- df_names %>%
# #     .[grepl("1", .)]
# #   #
# #   #
# #   # swap 1 and 2 in mod_1_names
# #   mod_2_names <- mod_1_names %>%
# #     stringr::str_replace_all("2", "3") %>%
# #     stringr::str_replace_all("1", "2") %>%
# #     stringr::str_replace_all("3", "1")
# #
# #   # mod_1_names <- df_names %>%
# #   #   .[grepl("samp_1", .)]
# #   #
# #   #
# #   # # swap 1 and 2 in mod_1_names
# #   # mod_2_names <- mod_1_names %>%
# #   #   stringr::str_replace_all("samp_1", "samp_2")
# #
# #
# #
# #   mod_1_names <- df_names %>%
# #     .[grepl("den_1", .)]
# #
# #
# #   # swap 1 and 2 in mod_1_names
# #   mod_2_names <- mod_1_names %>%
# #     stringr::str_replace_all("den_1", "den_2")
# #
# #   # with probability 0.5, swap samp_1 and samp_2
# #   # in randomly chosen rows of pp_df
# #   browser()
# #   nn <- nrow(pp_df)
# #   jj <- unique(sample(seq_len(nn), round(nn/2), replace = TRUE))
# #   tmp <- unname(pp_df[jj, mod_1_names])
# #   pp_df[jj, mod_1_names] <- unname(pp_df[jj, mod_2_names])
# #   pp_df[jj, mod_2_names] <- tmp
# #
# #   pp_df
# # }
# #
# #
# .swap_12 <- function(x) {
#   `%>%` <- magrittr::`%>%`
#   x %>%
#     stringr::str_replace_all("1", "3") %>%
#     stringr::str_replace_all("2", "1") %>%
#     stringr::str_replace_all("3", "2")
# }
#
# .rand_swap_mod_labs <- function(pp_df) {
#   `%>%` <- magrittr::`%>%`
#
#   pp_df1 <- data.table::copy(
#     pp_df
#   )[,
#     mod := ifelse(runif(nrow(pp_df)) <= 0.5, "mod1", "mod2")
#     ][,
#       den_curr := ifelse(mod == "mod1", den_mod1, den_mod2)
#       ]
#
#   pp_df1
# }
#
#
#
# .table_pp_sample_density2 <- function(n = 1000,
#                                       parlist1,
#                                       parlist2) {
#   allparlists <- list(`mod1` = parlist1,
#                       `mod2` = parlist2)
#
#   pp_samps <- purrr::map2(
#     allparlists,
#     names(allparlists),
#     function(parlist, ii) {
#       tmp <- do.call(
#         cbind,
#         .post_pred_sample_parlist(n, parlist)
#       )
#       tmp_dt <- as.data.frame(tmp)
#       tmp_dt$mod <- ii
#       tmp_dt
#     }
#   ) %>%
#     do.call(rbind, .)
#
#
#
#   allden_allmod <- lapply(
#     allparlists,
#     function(parlist) {
#       .eval_density_parlist(
#         parlist = parlist,
#         coords = data.matrix(pp_samps[, .(line, angle)]),
#         angle.max = 2*pi
#       )$denmat
#     }
#   )
#
#   names(allden_allmod) <- paste0("den_", names(allden_allmod))
#
#
#   pp_df <- cbind(pp_samps, as.data.frame(allden_allmod))
#
#   pp_df[,
#         `:=`(den_curr = ifelse(mod == "mod1", den_mod1, den_mod2),
#              den_half = (den_mod1 + den_mod2)/2)
#         ]
#
#   pp_df
#
# }
#
#
#
# .JS_pp_df2 <- function(pp_df) {
#
#   KL_mod_to_half <- pp_df[,
#                           log_den_diff :=
#                             log(den_curr) - log(den_half)
#                           ][,
#                             .(KL = mean(log_den_diff))
#                             ]
#   # [,
#   #   .(KL = mean(log_den_diff),
#   #     KL_var = var(log_den_diff)),
#   #   by = mod
#   #   ]
#
#   # N <- nrow(pp_df[mod == "mod1"])
#
#   JS <- JS_val <- sum(pmax(KL_mod_to_half$KL, 0))/2
#   # JS_se <- sqrt(sum(KL_mod_to_half$KL_var))/sqrt(N)
#   # attr(JS, "se") <- JS_se
#   # attr(JS, "lower") <- JS_val - 1.96 * JS_se
#   # attr(JS, "upper") <- JS_val + 1.96 * JS_se
#   #
#   #
#   # JSdist <- JSdist_val <- sqrt(JS_val/log(2))
#   # JSdist_se <- JS_se/(2 * sqrt(JS_val) * log(2))
#   # attr(JSdist, "se") <- JSdist_se
#   # attr(JSdist, "lower") <- JSdist_val - 1.96 * JSdist_se
#   # attr(JSdist, "upper") <- JSdist_val + 1.96 * JSdist_se
#   #
#   #
#   # list(divergence = JS,
#   #      distance = JSdist)
#   sqrt(JS/log(2))
#
# }
#
#
#
# # under equality of the two models, two random mixtures of the
# # predicitive distributions can have larger distance than the original
# # distributions. If the two original distributions are far from
# # each other, then it is unlikely that the two mixtures will have
# # a larger distance.
#
# # USING IMPORTANCE SAMPLING
# .JS_pp_df_randmix_importance <- function(pp_df) {
#
#   `:=` <- data.table::`:=`
#   pp_df1 <- data.table::copy(pp_df)
#
#   # randomly generate 2 mixing probabilities of mod1 and mod2
#   mix_prob <- runif(2)
#
#   # create two mixing densities -- these will play the role of
#   # new den1 and den2
#
#   # mix_den_i is evaluated at samples from mod_i
#   # This must be accounted for via importance
#   # sampling while computing Monte Carlo averaging
#
#
#   pp_df1[,
#          `:=`(
#            den_mix1 = mix_prob[1] * den_mod1 +
#              (1 - mix_prob[1]) * den_mod2,
#            den_mix2 = mix_prob[2] * den_mod1 +
#              (1 - mix_prob[2]) * den_mod2
#
#          )
#          ][,
#            `:=`(
#              den_mix_curr = ifelse(mod == "mod1",
#                                    den_mix1,
#                                    den_mix2),
#              den_mix_half = (den_mix1 + den_mix2)/2
#            )
#            ][,
#              imp_ratio := den_mix_curr/den_curr
#              ]
#
#   psis_imp_ratio <- lapply(
#     c("mod1", "mod2"),
#     function(mod_curr) {
#       suppressWarnings(
#         loo::psis(log(pp_df1[mod == mod_curr]$imp_ratio),
#                   r_eff = NA)$log_weights
#       )
#     }
#   ) %>%
#     unlist()
#
#   pp_df1$psis <- exp(psis_imp_ratio)
#
#
#
#   # need to multiply by importance sampling
#   # factor den_mix_curr/den_curr
#   KL_mod_to_half <- pp_df1[,
#                            .(KL = sum((log(den_mix_curr) -
#                                          log(den_mix_half)) *
#                                         psis) /
#                                sum(psis)
#
#                            ),
#                            by = mod
#                            ]
#
#
#   # if (any( KL_mod_to_half$KL < 0)) browser()
#
#   JS <- sum(pmax(KL_mod_to_half$KL, 0))/2
#
#   sqrt(JS/log(2))
# }
# #
# #
# #
# # # USING SAMPLES FROM TWO MIXTURES
# # .JS_pp_df_randmix <- function(pp_df) {
# #
# #   `:=` <- data.table::`:=`
# #   pp_df1 <- data.table::copy(pp_df)
# #
# #   # randomly generate 2 mixing probabilities of mod1 and mod2
# #   mix_prob <- runif(2)
# #
# #   browser()
# #
# #   # total number of samples from each model
# #   N <- nrow(pp_df1[mod == "mod1"])
# #
# #   is_mod1_mix1 <- (runif(N) <= mix_prob[1])
# #   is_mod1_mix2 <- (runif(N) <= mix_prob[2])
# #
# #   data.table::setkey(pp_df1, "mod") # arrange rows by mod1 then mod2
# #
# #   pp_df1[,
# #          `:=`(
# #            mix1_samp_indic = c(is_mod1_mix1, !is_mod1_mix1),
# #            mix2_samp_indic = c(is_mod1_mix2, !is_mod1_mix2)
# #          )
# #          ]
# #
# #
# #   KL <- rbind(pp_df1[mix1_samp_indic == TRUE,
# #                      .(den_mix1 = mix_prob[1] * den_mod1 +
# #                          (1 - mix_prob[1]) * den_mod2,
# #                        den_mix2 = mix_prob[2] * den_mod1 +
# #                          (1 - mix_prob[2]) * den_mod2,
# #                        mix_mod = "mix_mod1")
# #                      ],
# #               pp_df1[mix2_samp_indic == TRUE,
# #                      .(den_mix1 = mix_prob[1] * den_mod1 +
# #                          (1 - mix_prob[1]) * den_mod2,
# #                        den_mix2 = mix_prob[2] * den_mod1 +
# #                          (1 - mix_prob[2]) * den_mod2,
# #                        mix_mod = "mix_mod2")
# #                      ]
# #   )[,
# #     `:=`(
# #       den_mix_curr = ifelse(mix_mod == "mix_mod1",
# #                             den_mix1,
# #                             den_mix2),
# #       den_mix_half = (den_mix1 + den_mix2)/2
# #     )
# #     ][,
# #       .(KL = mean(log(den_mix_curr) -
# #                     log(den_mix_half))),
# #       by = mix_mod
# #       ]
# #
# #
# #   JS <- sum(pmax(KL$KL, 0))
# #
# #   sqrt(JS/log(2))
# # }
#
#  %>% JS_dist2 <- function(mod1, mod2,
#                           nsim = 1000,
#                           nperm = 1000) {
#
#   pp_df <- .table_pp_sample_density2(
#     n = nsim,
#     parlist1 = mod1,
#     parlist2 = mod2
#   )
#
#   obs <- .JS_pp_df2(pp_df)
#
#   # rand <- vapply(
#   #   seq_len(nperm),
#   #   function(ii) {
#   #     .JS_pp_df_randmix_importance(pp_df)
#   #   },
#   #   0
#   # )
#
#   # list(obs = obs,
#   #      rand = rand,
#   #      `P(rand > obs)` = mean(rand >= obs)
#   # )
#   obs
#
# }
#
