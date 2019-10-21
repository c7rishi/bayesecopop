#' Evaluate JS distance between two fitted (Dirichlet) mixture models
#' @param mod1,mod2 cyl_mix objects (outputs of fit_cyl_mix) whose distance is of interest.
#' @export
calc_JS_dist <- function(mod1, mod2,
                         nsim = 100,
                         nperm = 100,
                         show_progress_perm = FALSE) {
  stopifnot(is.cylmix(mod1), is.cylmix(mod2))

  parlist1 <- rstan::extract(mod1$stanmodel)
  parlist2 <- rstan::extract(mod2$stanmodel)

  cat("\n Computing observed distance..")
  pp_df <- .table_pp_sample_density(
    n = nsim,
    parlist1 = parlist1,
    parlist2 = parlist2
  )

  mod1_dens_cols <- names(pp_df) %>% .[grepl("den_indiv_mod1", .)]
  mod2_dens_cols <- names(pp_df) %>% .[grepl("den_indiv_mod2", .)]
  n_mod1_dens_cols <- length(mod1_dens_cols)
  n_mod2_dens_cols <- length(mod2_dens_cols)


  pp_df_list <- as.list(pp_df)[c(mod1_dens_cols, mod2_dens_cols)]

  full_col_sum <- Reduce(
    "+",
    pp_df_list
  )


  obs <- .JS_pp_df(
    pp_df_list = pp_df_list,
    mod1_dens_cols = 1:n_mod1_dens_cols,
    n_mod1_dens_cols = n_mod1_dens_cols,
    n_mod2_dens_cols = n_mod2_dens_cols,
    full_col_sum = full_col_sum,
    mod_names_vec = pp_df$mod,
    rand_swap = FALSE
  )


  cat("\nComputing distances under random permutations..")
  if (show_progress_perm) {
    pb <- txtProgressBar(0, nperm, style = 3)
  }


  rand <- vapply(
    seq_len(nperm),
    function(ii) {
      tmp <- .JS_pp_df(
        pp_df_list = pp_df_list,
        mod1_dens_cols = 1:n_mod1_dens_cols,
        n_mod1_dens_cols = n_mod1_dens_cols,
        n_mod2_dens_cols = n_mod2_dens_cols,
        full_col_sum = full_col_sum,
        mod_names_vec = pp_df$mod,
        rand_swap = TRUE
      )
      if (show_progress_perm) {
        setTxtProgressBar(pb, ii)
      }
      tmp
    },
    0
  )

  list(obs = obs,
       rand = rand,
       `P(rand >= obs)` = mean(rand >= obs)
  )
}

