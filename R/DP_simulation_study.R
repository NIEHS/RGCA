# RE_DP_simulation_study
get_emp_llh <- function(obs_res, boot_vals) {
  bw <- bw.SJ(boot_vals)
  # sum equivalent to a kernel density estimate at the point
  p_density <- sum(dnorm(obs_res, mean = boot_vals, sd = bw)) /
    length(boot_vals)
  return(log(p_density))
}

# a script to run a simulation where a method is chosen as the truth and the
# error from using a different method to make a prediction is computed
run_comparison_simulation <- function(save_plot = FALSE) {
  model_set <- c(
    "GCA (slope=1)",
    "RGCA (slope~1)",
    "RGCA (slope>>1)",
    "RGCA (2-step)",
    "RGCA (2-step KM)",
    "RGCA (sill<0)"
  ) # , "Synergy")
  models_to_run <- 1:6
  show_plots <- FALSE
  n_top <- 5
  mix_size <- c(10, 20)
  n_big_iters <- 50
  # record metrics
  record_MSE <- matrix(0,
    nrow = n_big_iters * length(models_to_run) *
      length(mix_size),
    ncol = 8
  )
  record_CRPS <- record_MSE
  record_MPD <- record_MSE

  set.seed(100)
  curr_idx <- 0
  for (mix_iter in 1:length(mix_size)) {
    for (outer_loop in 1:length(models_to_run)) {
      true_mix_model <- model_set[models_to_run[outer_loop]]

      n_chems <- mix_size[mix_iter]
      for (big_loop_iter in 1:n_big_iters) {
        # Generate individual curves:
        a_vec <- runif(n_chems, min = 1.5, max = 10)
        slope_vec <- runif(n_chems, min = .5, max = 1.5)
        ec_vec <- runif(n_chems, min = .1, max = 20)

        if (true_mix_model == "RGCA (sill<0)") {
          a_vec <- runif(n_chems, min = -10, max = 10)
        }
        if (true_mix_model == "RGCA (slope>>1)") {
          slope_vec <- runif(n_chems, min = .1, max = 10)
        }

        # generate mix concentrations:  each row of matrix is one mix dose
        # equipotent, equimolar, dominant
        n_samps <- 30
        equipot_conc_matrix <- matrix(0, nrow = n_samps, ncol = n_chems)
        equimol_conc_matrix <- matrix(0, nrow = n_samps, ncol = n_chems)

        for (chem_idx in 1:n_chems) {
          equipot_conc_matrix[, chem_idx] <-
            ec_vec[chem_idx] / (10^seq(3, -2, length.out = n_samps))
          equimol_conc_matrix[, chem_idx] <-
            1 / (10^seq(4, -2, length.out = n_samps))
        }

        # Default case:  RGCA is true model
        sim_param_matrix <- as.matrix(cbind(
          "a" = a_vec,
          "b" = ec_vec,
          "c" = slope_vec,
          "max_R" = max(a_vec),
          "d" = 0
        ))
        mix_calc <- mix_function_generator(
          sim_param_matrix,
          rep(1, n_chems)
        )

        if (true_mix_model == "GCA (slope=1)") {
          slope_vec <- slope_vec * 0 + 1
          sim_param_matrix <- as.matrix(cbind(
            "a" = a_vec,
            "b" = ec_vec,
            "c" = slope_vec,
            "max_R" = max(a_vec),
            "d" = 0
          ))
          mix_calc <- mix_function_generator(
            sim_param_matrix,
            rep(1, n_chems)
          )
        }
        if (true_mix_model == "Synergy") {
          sim_param_matrix <- as.matrix(cbind(
            "a" = a_vec,
            "b" = ec_vec,
            "c" = slope_vec,
            "max_R" = max(a_vec),
            "d" = 0
          ))
          mix_calc <- mix_function_generator(sim_param_matrix,
            rep(1, n_chems),
            synergy_const = -1000
          )
        }
        # Case 3: RGCA with two step model
        if (true_mix_model == "RGCA (2-step)") {
          # assume either random clustering or use DP
          clust_rand <- random_assignment(n_chems, sample(2:3, 1))
          # generate mixture curve from IA/GCA
          sim_param_matrix <- as.matrix(cbind(
            "a" = a_vec,
            "b" = ec_vec,
            "c" = slope_vec,
            "max_R" = max(a_vec),
            "d" = 0
          ))
          mix_calc <- mix_function_generator(
            sim_param_matrix,
            clust_rand
          )
        }

        kmeans_mat <- as.matrix(cbind(
          "a" = a_vec,
          "b" = ec_vec,
          "c" = slope_vec
        ))
        KM_clust <- kmeans(kmeans_mat, centers = 3)$cluster
        if (true_mix_model == "RGCA (2-step KM)") {
          mix_calc <- mix_function_generator(sim_param_matrix, KM_clust)
        }


        # compute true curves given mix calc
        true_response_equipot <-
          sapply(1:n_samps, function(x) mix_calc(equipot_conc_matrix[x, ]))

        if (show_plots) {
          x_range <- range(c(
            rowSums(equimol_conc_matrix),
            rowSums(equipot_conc_matrix)
          ))
          plot(rowSums(equipot_conc_matrix),
            true_response_equipot,
            log = "x",
            ylim = c(0, max(a_vec)),
            xlim = x_range
          )
          true_response_equimol <-
            sapply(1:n_samps, function(x) mix_calc(equimol_conc_matrix[x, ]))
          points(rowSums(equimol_conc_matrix), true_response_equimol, col = 2)
        }

        # compare GCA, IA, our method via UQ: CRPS, etc
        sim_fixed_params <- list(
          "sill_params" = a_vec,
          "sill_sd" = a_vec * 0,
          "ec50_params" = ec_vec,
          "ec50_stdev" = ec_vec * 0,
          "u_RE_sd_params" = rep(0, n_chems),
          "v_RE_sd_params" = rep(0, n_chems),
          "slope_params" = slope_vec
        )


        # GCA and CA
        GCA_assign <- 1
        names(GCA_assign) <- do.call(paste, as.list(rep(1, n_chems)))
        GCA_assign_vec <-
          as.numeric(strsplit(names(GCA_assign), split = " ")[[1]])
        param_matrix_GCA <- as.matrix(cbind(
          "a" = a_vec,
          "b" = ec_vec,
          "c" = slope_vec * 0 + 1,
          "max_R" = max(a_vec)
        ))
        sim_mix_funs_GCA <- mix_function_generator(
          param_matrix_GCA,
          GCA_assign_vec
        )
        sim_mix_funs_CA <- mix_function_generator(param_matrix_GCA,
          GCA_assign_vec,
          scale_CA = TRUE
        )

        # RGCA
        param_matrix_RGCA <- as.matrix(cbind(
          "a" = a_vec,
          "b" = ec_vec,
          "c" = slope_vec,
          "max_R" = max(a_vec)
        ))
        sim_mix_funs_RGCA <- mix_function_generator(
          param_matrix_RGCA,
          GCA_assign_vec
        )
        # IA
        IA_assign <- 1
        names(IA_assign) <- do.call(paste, as.list(1:n_chems))
        IA_assign_vec <-
          as.numeric(strsplit(names(IA_assign), split = " ")[[1]])
        sim_mix_funs_IA <-
          mix_function_generator(param_matrix_RGCA, IA_assign_vec)

        # RGCA 2-STEP, random clusters
        n_bootstraps <- 10
        # generate some random cluster assignments
        rand_clust_mat <- rbind(
          rep(1, n_chems),
          t(sapply(1:4, FUN = function(x) {
            random_assignment(n_chems, 2 + 0 * x)
          })),
          t(sapply(1:4, FUN = function(x) {
            random_assignment(n_chems, 3 + 0 * x)
          })),
          1:n_chems
        )

        rand_clust_assign <- rep(1, nrow(rand_clust_mat))
        names(rand_clust_assign) <-
          apply(rand_clust_mat,
            MARGIN = 1,
            FUN = function(rx) do.call(paste, as.list(rx))
          )
        randclust_par_list <- list(
          "cluster_assign" = rand_clust_assign,
          "centers" = slope_vec,
          "cent_sd" = rep(0, n_chems)
        )
        RGCA_randclust_par_list <- c(sim_fixed_params, randclust_par_list)

        joint_model_fun_rand <- function(x) {
          create_mix_calc_from_summary(x, RGCA_randclust_par_list)
        }
        sim_mix_funs_2step <- sapply(1:n_bootstraps,
          FUN = function(x) joint_model_fun_rand(x)
        )

        if (FALSE) {
          # RGCA 2-STEP, DP on sill
          cluster_prob <- clust_centers_w_prob$probs
          cluster_assign <- clust_centers_w_prob$assign
          samp_idx <- sample(1:n_top,
            size = n_bootstraps,
            prob = cluster_prob, replace = TRUE
          )
          cluster_par_list_sim_sill <- c(
            sim_fixed_params,
            list(
              "cluster_assign" = cluster_assign,
              "centers" = slope_vec,
              "cent_sd" = rep(0, n_chems)
            )
          )
          joint_model_fun_sill <- function(x) {
            create_mix_calc_from_summary(x, cluster_par_list_sim_sill)
          }
          sim_mix_funs_2step_DP <-
            sapply(samp_idx, FUN = function(x) joint_model_fun_sill(x))
        }

        sim_mix_funs_2step_KM <- mix_function_generator(
          param_matrix_RGCA,
          KM_clust
        )
        ## Plot results for simulated data
        conc_mat <- equipot_conc_matrix
        curve_data_2step <- matrix(0, ncol = n_samps, nrow = n_bootstraps)
        curve_data_2stepKM <- matrix(0, ncol = n_samps, nrow = 1)
        curve_data_RGCA <- matrix(0, ncol = n_samps, nrow = 1)
        curve_data_GCA <- matrix(0, ncol = n_samps, nrow = 1)
        curve_data_CA <- matrix(0, ncol = n_samps, nrow = 1)
        curve_data_IA <- matrix(0, ncol = n_samps, nrow = 1)
        for (row_idx in 1:n_samps) {
          conc_val <- conc_mat[row_idx, ] #+1e-40
          curve_data_2step[, row_idx] <- sapply(sim_mix_funs_2step,
            FUN = function(x) x(conc_val)
          )
          curve_data_2stepKM[, row_idx] <- sim_mix_funs_2step_KM(conc_val)
          curve_data_RGCA[, row_idx] <- sim_mix_funs_RGCA(conc_val)
          curve_data_GCA[, row_idx] <- sim_mix_funs_GCA(conc_val)
          curve_data_CA[, row_idx] <- sim_mix_funs_CA(conc_val)
          curve_data_IA[, row_idx] <- sim_mix_funs_IA(conc_val)
        }

        Cx_axis_values <- rowSums(conc_mat)

        if (show_plots) {
          make_ggplot(
            Cx_axis_values,
            true_response_equipot,
            curve_data_RGCA,
            curve_data_GCA,
            t(curve_data_2step[3, ])
          )
        }

        # get MSE, Logscore, CRPS  ####
        n_x_values <- ncol(Cx)
        RGCA_llh <- 0
        RGCA_mse <- 0
        RGCA_crps <- 0
        GCA_llh <- 0
        GCA_mse <- 0
        GCA_crps <- 0
        IA_llh <- 0
        IA_mse <- 0
        IA_crps <- 0
        RGCA_2step_llh <- 0
        RGCA_2step_mse <- 0
        RGCA_2step_crps <- 0
        RGCA_2stepKM_llh <- 0
        RGCA_2stepKM_mse <- 0
        RGCA_2stepKM_crps <- 0
        CA_mse <- 0
        CA_crps <- 0

        # average % difference from EC50?
        RGCA_mpd <- 0
        GCA_mpd <- 0
        CA_mpd <- 0
        IA_mpd <- 0
        RGCA_2step_mpd <- 0
        RGCA_2stepKM_mpd <- 0

        # loop through each dose and get pointwise scores across all samples
        for (x_idx in 1:n_x_values) {
          ClGCA_llh <- ClGCA_llh + get_emp_llh(
            true_response_equipot[x_idx],
            curve_data[, x_idx]
          )
          GCA_llh <- GCA_llh + get_emp_llh(
            true_response_equipot[x_idx],
            curve_data_GCA[, x_idx]
          )
          IA_llh <- IA_llh + get_emp_llh(
            true_response_equipot[x_idx],
            curve_data_IA[, x_idx]
          )
          RGCA_2step_crps <- RGCA_2step_crps +
            scoringRules::crps_sample(
              true_response_equipot[x_idx],
              curve_data_2step[, x_idx]
            )
          RGCA_2stepKM_crps <- RGCA_2stepKM_crps +
            scoringRules::crps_sample(
              true_response_equipot[x_idx],
              curve_data_2stepKM[, x_idx]
            )
          RGCA_crps <- RGCA_crps +
            scoringRules::crps_sample(
              true_response_equipot[x_idx],
              curve_data_RGCA[, x_idx]
            )
          GCA_crps <- GCA_crps +
            scoringRules::crps_sample(
              true_response_equipot[x_idx],
              curve_data_GCA[, x_idx]
            )
          CA_crps <- CA_crps +
            scoringRules::crps_sample(
              true_response_equipot[x_idx],
              curve_data_CA[, x_idx]
            )
          IA_crps <- IA_crps +
            scoringRules::crps_sample(
              true_response_equipot[x_idx],
              curve_data_IA[, x_idx]
            )

          RGCA_2step_mse <- RGCA_2step_mse +
            (true_response_equipot[x_idx] -
              median(curve_data_2step[, x_idx]))^2 / n_x_values
          RGCA_2stepKM_mse <- RGCA_2stepKM_mse +
            (true_response_equipot[x_idx] -
              median(curve_data_2stepKM[, x_idx]))^2 / n_x_values
          RGCA_mse <- RGCA_mse +
            (true_response_equipot[x_idx] -
              median(curve_data_RGCA[, x_idx]))^2 / n_x_values
          GCA_mse <- GCA_mse +
            (true_response_equipot[x_idx] -
              median(curve_data_GCA[, x_idx]))^2 / n_x_values
          CA_mse <- CA_mse + (true_response_equipot[x_idx] -
            median(curve_data_CA[, x_idx]))^2 / n_x_values
          IA_mse <- IA_mse + (true_response_equipot[x_idx] -
            median(curve_data_IA[, x_idx]))^2 / n_x_values

          RGCA_2step_mpd <- RGCA_2step_mpd +
            (true_response_equipot[x_idx] -
              median(curve_data_2step[, x_idx])) / true_response_equipot[x_idx]
          RGCA_2stepKM_mpd <- RGCA_2stepKM_mpd +
            (true_response_equipot[x_idx] -
              median(curve_data_2stepKM[, x_idx])) /
              true_response_equipot[x_idx]
          RGCA_mpd <- RGCA_mpd +
            (true_response_equipot[x_idx] -
              median(curve_data_RGCA[, x_idx])) / true_response_equipot[x_idx]
          GCA_mpd <- GCA_mpd +
            (true_response_equipot[x_idx] -
              median(curve_data_GCA[, x_idx])) / true_response_equipot[x_idx]
          CA_mpd <- CA_mpd +
            (true_response_equipot[x_idx] -
              median(curve_data_CA[, x_idx])) / true_response_equipot[x_idx]
          IA_mpd <- IA_mpd +
            (true_response_equipot[x_idx] -
              median(curve_data_IA[, x_idx])) / true_response_equipot[x_idx]
        }

        curr_idx <- curr_idx + 1
        record_MSE[curr_idx, ] <- c(
          RGCA_2step_mse,
          RGCA_2stepKM_mse,
          RGCA_mse,
          GCA_mse,
          CA_mse,
          IA_mse,
          true_mix_model,
          n_chems
        )

        record_MPD[curr_idx, ] <- c(
          RGCA_2step_mpd,
          RGCA_2stepKM_mpd,
          RGCA_mpd,
          GCA_mpd,
          CA_mpd,
          IA_mpd,
          true_mix_model,
          n_chems
        )

        record_CRPS[curr_idx, ] <- c(
          RGCA_2step_crps,
          RGCA_2stepKM_crps,
          RGCA_crps,
          GCA_crps,
          CA_crps,
          IA_crps,
          true_mix_model,
          n_chems
        )

        # annotation offset: 10% of max(curve_data_IA)
      }
    }
  }
  beepr::beep()

  if (save_plot) {
    pdf(file = "sim_study_invertx.pdf", width = 7, height = 5)
    data_col_names <- c("RGCA (2-step)", "RGCA (2-step KM)", "RGCA", "GCA", "CA", "IA", "Mix ID", "Num Chems")
    sim_study_boxplot(record_MSE, data_col_names)
    dev.off()
  }
}
