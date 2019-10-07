library(bayesecopop)
# create dataset
set.seed(100)
dat_new <- gen_dat(ndata = 500, ncomp = 4)
dat_new$pars[-6]

plot(dat_new$dat)


# plot the contour at a specific parameter value
bayesecopop:::contour_plot_singlepar(dat_new$dat,
                                     parlist = dat_new$pars,
                                     alpha = 0.3)

dd <- dat_new$dat


# fit_vb_finite <- fit_cyl_mix(
#   dd,
#   ncomp = 4,
#   line_col = "time",
#   angle_col = "angle",
#   sampling_method = "vb",
#   output_samples = 5000
# )

fit_vb_dirichlet <- fit_cyl_mix(
  dd,
  mixture_type = "dirichlet",
  ncomp = 6,
  line_col = "time",
  angle_col = "angle",
  sampling_method = "vb",
  output_samples = 1000,
)
estimate_nonempty_ncomp(fit_vb_dirichlet)


# fit_vb_finite_sparse <- fit_cyl_mix(
#   dd,
#   mixture_type = "finite",
#   overfit = "leave-empty",
#   line_col = "time",
#   angle_col = "angle",
#   alpha0_finite =  1e-6,
#   sampling_method = "vb",
#   output_samples = 5000
# )
# fit_mcmc_finite_sparse <- fit_cyl_mix(
#   dd,
#   mixture_type = "finite",
#   overfit = "leave-empty",
#   line_col = "time",
#   angle_col = "angle",
#   alpha0_finite =  1e-6,
#   sampling_method = "mcmc"
# )

fit_mcmc_dirichlet <- fit_cyl_mix(
  dd,
  mixture_type = "dirichlet",
  ncomp = 6,
  line_col = "time",
  angle_col = "angle",
  sampling_method = "mcmc"
)



estimate_nonempty_ncomp(fit_mcmc_finite_sparse)



# fits <- lapply(2:8,
#                function(j)
#                {
#                  set.seed(100)
#                  cat("\n")
#                  cat("**************************************\n")
#                  cat("ncomp =", j, "\n")
#                  cat("**************************************\n")
#                  fit_mixmod_stanvb(dat_new$dat, ncomp = j)
#
#                })
#
#
#
# tmp <- fit_incremental_mixmod(dat_new$dat)
#
# fits_stan <- lapply(fits, "[[", "st")
# fits_map <- lapply(fits, "[[", "init")

# names(fits_vb) <- names(fits_map) <-
#   paste("m", 2:8, sep = " = ")


mm <- post_pred_density(fit_mcmc, data = dat_new$dat,
                        line_col = "time")


# plot contour by average density from all posterior samples
bayesecopop:::contour_plot_stan(
  fit_vb_finite$stanmodel,
  data = dat_new$dat,
  alpha = 0.3)


bayesecopop:::contour_plot_stan(
  fit_vb_dirichlet$stanmodel,
  data = dat_new$dat,
  alpha = 0.3)

bayesecopop:::contour_plot_stan(
  fit_vb_finite_sparse$stanmodel,
  data = dat_new$dat,
  alpha = 0.3)

loo_fits <- lapply(fits_vb, loo)

loo::compare(x = loo_fits)



# overfit with leave some group empty
set.seed(100)
fit_10 <- fit_mixmod_stanvb(dat_new$dat, ncomp = 30,
                            overfit = "leave-empty")
contour_plot_singlepar(dat_new$dat, parlist = fit_10$init)

contour_plot_stan(fit_10$stan_vb,
                  dat_new$dat)

est_ncomp_stan(fit_10$stan_vb)
