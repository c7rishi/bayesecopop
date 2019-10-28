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





# computes Jensen-Shannon distance
# (in [0, 1]) scale, with permutation
# p-values
.JS_pp_df <- function(pp_df_list,
                      mod1_dens_cols,
                      n_mod1_dens_cols,
                      n_mod2_dens_cols,
                      full_col_sum,
                      mod_names_vec,
                      rand_swap = FALSE) {

  if (rand_swap) {
    # browser()
    # mix_prop <- runif(2)
    # pp_df_list1 <- pp_df_list
    # pp_df_list <- lapply(
    #   mix_prop,
    #   function(this_mix_prop) {
    #     mapply(
    #       function(x, y) {
    #         this_mix_prop * x + (1-this_mix_prop) * y
    #       },
    #       pp_df_list[mod1_dens_cols],
    #       pp_df_list[-mod1_dens_cols],
    #       SIMPLIFY = FALSE
    #     )
    #   }
    # ) %>%
    #   unlist(recursive = FALSE)

    # mod_names_vec <- sample(mod_names_vec)

    #
    #     mod1_dens_cols <- sample(1:(n_mod1_dens_cols+n_mod1_dens_cols),
    #                              n_mod1_dens_cols, replace = TRUE)
    #     mod2_dens_cols <- sample(1:(n_mod1_dens_cols+n_mod1_dens_cols),
    #                              n_mod1_dens_cols, replace = TRUE)

    # mod_names_vec <- sample(mod_names_vec)

    # label each density evaluation to be from
    # model 1 with probability mix_prop
    # mod1_dens_cols <- which(
    #   runif(n_mod1_dens_cols+n_mod1_dens_cols) <= mix_prop
    # )
    # mod_names_vec <- c("mod1", "mod2")[
    #   1 + (runif(length(mod_names_vec)) <= mix_prop)
    #   ]
    # mod1_dens_cols <- sample(1:(n_mod1_dens_cols+n_mod1_dens_cols),
    #                          n_mod1_dens_cols)
    # mod_names_vec <- sample(mod_names_vec)

    # browser()
    pp_rows <- lapply(
      unique(mod_names_vec),
      function(this_mod_name) {
        sample(
          seq_len(length(mod_names_vec)),
          sum(mod_names_vec == this_mod_name),
          replace = TRUE
        )
      }
    ) %>%
      unlist()

    pp_df_list <- lapply(
      pp_df_list,
      "[",
      pp_rows
    )


  }


  tmp_sum <- Reduce(
    "+",
    pp_df_list[mod1_dens_cols]
  )


  den_mod1 <- tmp_sum/n_mod1_dens_cols
  # den_mod2 <- (full_col_sum - tmp_sum)/n_mod2_dens_cols
  if (rand_swap) {
    den_mod2 <- Reduce(
      "+",
      pp_df_list[-mod1_dens_cols]
    ) / n_mod2_dens_cols
  } else {
    den_mod2 <- (full_col_sum - tmp_sum)/n_mod2_dens_cols
  }

  mod1_indic <- mod_names_vec == "mod1"
  den_curr <- ifelse(mod1_indic,
                     den_mod1, den_mod2)
  den_half <- (den_mod1 + den_mod2)/2
  # if (!rand_swap) {
  #   den_half <- (den_mod1 + den_mod2)/2
  # } else {
  #   pmix <- runif(1)
  #   den_half <- pmix * den_mod1 + (1 - pmix) * den_mod2
  # }

  log_den_diff <- log(den_curr) - log(den_half)

  mod1_idx <- which(mod1_indic)
  KL_mod_to_half <- list(log_den_diff[mod1_idx],
                         log_den_diff[-mod1_idx]) %>%
    lapply(mean) %>%
    unlist()

  JS <- sum(pmax(KL_mod_to_half, 0))/2

  sqrt(JS/log(2))

}


