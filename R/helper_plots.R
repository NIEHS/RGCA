## Plotting functions for DP+RGCA method

library(ggplot2)
library(cowplot)
library(reshape2)


#' Plot individual dose response curves
#'
#' Plots a few hard-coded samples of the data for a comparison used in the
#' manuscript.  Shows why a random effect model is needed for the Tox21 data.
#' Requires the data cleaning script to be previously sourced,
#' source("Code/tox21_prep_data.R").
#'
#' @param df a data frame with the prepared tox21 data from
#'   source("Code/tox21_prep_data.R")
#'
#' @return a ggplot with four examples of the tox21 data
#' @export
#'
#' @examples
plot_sample_indv <- function(df) {
  chems_to_plot <- c("Hydroxyflutamide",
                     "Progesterone",
                     "Fluoxymestrone",
                     "Ethylenediamine")
  # get the names of all chemicals in the data
  sample_names <- unique(pure_df$Sample.Name)
  # take indices from Global env, set during data prep
  conc_idx <- .GlobalEnv$conc_idx
  resp_idx <- .GlobalEnv$resp_idx
  # extract x-y values and a map of replicate to idx
  x_vals <- c()
  y_vals <- c()
  for (chem_nm in chems_to_plot) {
    cm_idx <- which(df$Sample.Name == chem_nm)
    # add CAS to chem label for clarity
    chem_CAS <- names(chem_map_plain)[which(chem_map_plain == chem_nm)]
    chem_label <- paste(chem_nm, " (", chem_CAS, ")", sep = "")
    # create long-format data frame with concentrations (x-values)
    new_x <- cbind(reshape2::melt(t(df[cm_idx, conc_idx])), "chem" = chem_label)
    # substitute replicate index with replicate number
    new_x$Var2 <- sapply(new_x$Var2, FUN = function(x) which(x == cm_idx))
    # create long-format data with responses (y-value)
    new_y <- cbind(reshape2::melt(t(df[cm_idx, resp_idx])))
    x_vals <- rbind(x_vals, new_x)
    y_vals <- rbind(y_vals, new_y)
  }
  # combine into a data frame and subset to drop duplicates
  gg_df <- cbind(x_vals, y_vals)[, c(4, 2, 3, 7)]
  names(gg_df) <- c("chem", "idx", "Concentration", "Response")
  ggplot2::ggplot(data = gg_df, aes(Concentration,
                                    Response,
                                    colour = as.factor(idx))) +
    geom_point() +
    geom_line() +
    scale_x_continuous(trans = "log10") +
    theme(legend.position = "top") +
    facet_wrap(vars(chem), scales = "free") +
    guides(colour = guide_legend(nrow = 1)) +
    labs(colour = "Tox21 AR-luc Replicates")
}
#nolint pdf("Output/replicate_ex.pdf", width =7, height = 4)
#nolint plot_sample_indv(df)
#nolint dev.off()



# Mixture DRC's ####

#' Makes a ggplot of the mixture response curves
#'
#' Creates plots showing the various methods as banded curves where the bands
#' are the posterior credible intervals and the observed mixture response as
#' points.
#'
#' @param Cx_axis_values matrix of mixture concentrations
#' @param resp_y_values matrix of observed mixture responses
#' @param curve_data matrix of predicted mixture responses, with a row for each
#'   bootstrapped calculator, RGCA method
#' @param curve_data_GCA matrix of predicted mixture responses, GCA method
#' @param curve_data_IA matrix of predicted mixture responses, IA method
#'
#' @return a ggplot
#' @export
make_ggplot <- function(Cx_axis_values,
                        resp_y_values,
                        curve_data,
                        curve_data_GCA = NULL,
                        curve_data_IA = NULL) {
  curve_data_quantiles <- t(apply(curve_data,
                                  MARGIN = 2,
                                  FUN = function(x) {
                                    quantile(x,
                                             probs = c(0.05, .5, .95),
                                             na.rm = TRUE)
                                  }))
  gg_df <- data.frame(Cx_axis_values, resp_y_values, curve_data_quantiles)
  names(gg_df) <- c("Cx", "resp_y", "low", "med", "up")
  hill_pred_plot <- ggplot2::ggplot(gg_df, aes(Cx, resp_y)) +
    geom_point(shape = 21,
               colour = "black",
               fill = "white",
               size = 3,
               stroke = 1.5) +
    scale_x_continuous(trans = "log10") +
    geom_line(aes(Cx, med, colour = "RE+DP")) +
    geom_ribbon(aes(ymin = low, ymax = up, fill = "RE+DP"),
                alpha = 0.1,
                color = "black", linetype = "dotted")
  if (!is.null(curve_data_GCA)) {
    curve_GCA_data_quantiles <-
      t(apply(curve_data_GCA,
              MARGIN = 2,
              FUN = function(x) {
                quantile(x, probs = c(0.05, .5, .95), na.rm = TRUE)
              }))
    GCA_df <- data.frame(Cx_axis_values,
                         resp_y_values,
                         curve_GCA_data_quantiles)
    names(GCA_df) <- c("Cx", "resp_y", "GCA_low", "GCA_med", "GCA_up")
    hill_pred_plot <- hill_pred_plot +
      geom_line(data = GCA_df, aes(Cx, GCA_med, colour = "GCA")) +
      geom_ribbon(
        data = GCA_df, aes(ymin = GCA_low, ymax = GCA_up, fill = "GCA"),
        alpha = 0.1,
        color = "black", linetype = "dotted"
      ) +
      scale_colour_manual(name = "Legend",
                          values = c("RE+DP" = "black", "GCA" = "red")) +
      scale_fill_manual(name = "Legend",
                        values = c("RE+DP" = "green", "GCA" = "red"))
  }
  if (!is.null(curve_data_IA)) {
    curve_IA_data_quantiles <- t(apply(curve_data_IA,
                                       MARGIN = 2,
                                       FUN = function(x) {
                                         quantile(x,
                                                  probs = c(0.05, .5, .95),
                                                  na.rm = TRUE)
                                       }))
    IA_df <- data.frame(Cx_axis_values, resp_y_values, curve_IA_data_quantiles)
    names(IA_df) <- c("Cx", "resp_y", "IA_low", "IA_med", "IA_up")
    hill_pred_plot <- hill_pred_plot +
      geom_line(data = IA_df, aes(Cx, IA_med, colour = "IA")) +
      geom_ribbon(
        data = IA_df, aes(ymin = IA_low, ymax = IA_up, fill = "IA"),
        alpha = 0.1,
        color = "black", linetype = "dotted"
      ) +
      scale_colour_manual(name = "Legend",
                          values = c("RE+DP" = "black",
                                     "GCA" = "red",
                                     "IA" = "blue")) +
      scale_fill_manual(name = "Legend",
                        values = c("RE+DP" = "green",
                                   "GCA" = "red",
                                   "IA" = "blue"))
  }
  hill_pred_plot
}


#' Makes a ggplot of the mixture response curves
#'
#' Creates plots showing the various methods as banded
#' curves where the bands are the posterior credible intervals and the observed
#' mixture response as points.  Slightly more general than the make_ggplot
#' function because a list of matrices can be passed rather than a separate
#' matrix for each prediction method.
#'
#' @param Cx_axis_values matrix of mixture concentrations
#' @param curve_data_list a list of matrices of predicted mixture responses.
#'   The rows of the matrices correspond to bootstrapped calculators, and each
#'   matrix can represent a different method.
#' @param resp_y_values matrix of observed mixture responses
#'
#' @return a ggplot
#' @export
#'
make_ggplot_many <- function(Cx_axis_values, resp_y_values, curve_data_list) {
  # plt_id specified from the name of the curve calculator
  plt_id <- names(curve_data_list)
  num_curves <- length(curve_data_list)
  if (num_curves == 0) {
    return(NA)
  }
  if (length(plt_id) != num_curves) print("curves not named")
  quantile_df_colnames <- c("Cx", "resp_y", "low", "med", "up", "ID")
  curve_data_quantiles <- t(apply(curve_data_list[[1]],
                                  MARGIN = 2,
                                  FUN = function(x) {
                                    quantile(x,
                                             probs = c(0.05, .5, .95),
                                             na.rm = TRUE)
                                  }))
  gg_df <- data.frame(Cx_axis_values,
                      resp_y_values,
                      curve_data_quantiles,
                      plt_id[1])
  names(gg_df) <- quantile_df_colnames
  gg_df$Replicate <- 1
  hill_pred_plot <- ggplot2::ggplot(gg_df, aes(Cx, resp_y)) +
    geom_point(shape = 21, colour = "black",
               fill = "white",
               size = 0,
               stroke = 1) +
    scale_x_continuous(trans = "log10") +
    geom_line(aes(Cx, med, colour = ID)) +
    geom_ribbon(aes(ymin = low,
                    ymax = up,
                    fill = plt_id[1]),
                alpha = 0.1, # if needed: fill = "green"
                color = "black",
                linetype = "dotted")
  for (curve_idx in 2:num_curves) {
    curve_data_quantiles <- t(apply(curve_data_list[[curve_idx]],
                                    MARGIN = 2,
                                    FUN = function(x) {
                                      quantile(x,
                                               probs = c(0.05, .5, .95),
                                               na.rm = TRUE)
                                    }))
    gg_df <- data.frame(Cx_axis_values,
                        resp_y_values,
                        curve_data_quantiles,
                        plt_id[curve_idx])
    names(gg_df) <- quantile_df_colnames
    hill_pred_plot <- hill_pred_plot +
      geom_line(data = gg_df, aes(Cx, med, colour = ID)) +
      geom_ribbon(data = gg_df,
                  aes(ymin = low, ymax = up, fill = ID),
                  alpha = 0.1,
                  color = "black",
                  linetype = "dotted")
  }
  # coloring:  assume first plot is our method, use black line with green curve
  line_colors <- c("black", "red", "blue", "orange", "magenta", "yellow")
  fill_colors <- c("green", "red", "blue", "orange", "magenta", "yellow")
  legend_color_values <- 1:num_curves
  fill_vals <- 1:num_curves
  fill_vals[1] <- 3
  if (num_curves > 2) {
    legend_color_values <- 1:(num_curves + 1)
    legend_color_values <- legend_color_values[-3]
    fill_vals <- 1:(num_curves + 1)
    fill_vals[1] <- 3
    fill_vals <- fill_vals[-3]
  }
  if (num_curves < length(line_colors)) {
    legend_color_values <- line_colors[1:num_curves]
    fill_vals <- fill_colors[1:num_curves]
  }
  names(legend_color_values) <- plt_id
  names(fill_vals) <- plt_id
  hill_pred_plot <- hill_pred_plot +
    scale_colour_manual(name = "Legend", values = legend_color_values) +
    scale_fill_manual(name = "Legend", values = fill_vals)
  return(hill_pred_plot)
}

#' Wrapper function to plot dummy concentrations
#'
#' A dummy concentration refers to a mixture with only one chemical.  This
#' plotting function is just used to test that the random effect model produces
#' a good fit.
#'
#' @param Cx matrix of concentrations or doses
#' @param y_i matrix of predicted responses
#' @param tot_par_list list of fitted parameters
#' @param replicate_sets a list of vectors.  The vectors contain the indices of
#'   replicates for a particular chemical.
#' @param bootstrapped_calc_list a list of lists of predictors with bootstrapped
#'   parameters.
#' @param test_idx the index of the individual chemical to plot
#'
#' @return
#' @export
#'
#' @examples
plot_dummy_mixture <- function(Cx,
                               y_i,
                               tot_par_list,
                               replicate_sets,
                               bootstrapped_calc_list,
                               test_idx = 1) {
  # test along original chemical concentrations (equal dose)
  n_dose <- ncol(Cx)
  dummy_conc_mat <- matrix(0, nrow = n_dose, ncol = n_chems)
  dummy_conc_mat[, test_idx] <- Cx[test_idx, ]
  curve_list <- predict_mix_response_many(n_dose = n_dose,
                                          chem_conc_matr = dummy_conc_mat,
                                          bootstrap_calc_list =
                                            bootstrapped_calc_list,
                                          default_entry = NA)
  specific_idx <- replicate_sets[[test_idx]][1]
  obs_resp <- y_i[specific_idx, ]
  make_ggplot_many(Cx[specific_idx, ], obs_resp, curve_list)
}

#' Plot the mixture response
#'
#' This method creates the manuscript plots that show the predicted median
#' mixture response, the credible interval bands around the median response, and
#' the observed mixture response with replicates.
#'
#' @param plot_set a vector of indices that specify which mixtures are predicted
#'   and plotted.  The indices correspond to the observed mixture response data
#'   file.
#' @param mix_df A data frame with the observed mixture response data
#' @param mix_conc_df a data frame with the concentrations of the mixture (not
#'   split into constituents)
#' @param mix_guide A data frame that specifies the concentrations of the
#'   constituent chemicals for each mixture at the highest dose, along with name
#'   and other meta data.
#' @param bootstrap_calc_list a list of functions that predict a mixture
#'   response given a concentration.  The functions are created by the factory
#'   function mix_function_generator
#'
#' @return
#' @export
#'
#' @examples
plot_mixture_response <- function(plot_set, mix_df, mix_conc_df, mix_guide,
                                  bootstrap_calc_list) {
  resp_idx <- .GlobalEnv$resp_idx
  n_calcs <- length(bootstrap_calc_list)
  # preallocate matrix for score record
  score_record <- matrix(0,
                         nrow = sum(mix_df[plot_set, ]$ReplicateSet == 1),
                         ncol = 1 + n_calcs * 3)
  for (iter_idx in seq(length(plot_set))) {
    mix_idx <- plot_set[iter_idx]
    chem_conc_matr <- get_conc_matrix(mix_idx)
    # first 15 dose values are real, second 15 are extrapolation.
    # Here, just extrapolating 10 less doses
    n_dose <- nrow(chem_conc_matr) - 10
    chem_conc_matr_reordered <- chem_conc_matr[, pure_ordering]
    curve_data_list <- predict_mix_response_many(n_dose,
                                                 chem_conc_matr_reordered,
                                                 bootstrap_calc_list)
    Cx_axis_values <- array(unlist(mix_conc_df[mix_idx, ]))[1:n_dose]
    mix_replicates <- which(mix_df$CAS == mix_df$CAS[mix_idx])
    resp_y_values <- array(unlist(mix_df[mix_idx, resp_idx]))
    # allow extrapolation; match dim of Cx by adding null values to y observed
    y_response_ext <- c(resp_y_values, rep(NA, n_dose - length(resp_y_values)))
    myplot <- make_ggplot_many(Cx_axis_values,
                               y_response_ext,
                               curve_data_list)
    # add replicates to single plot
    other_repl_idx <- which(mix_df$CAS == mix_df$CAS[mix_idx])
    unplotted_repl <- other_repl_idx
    y_unplt_repl <- rbind(t(mix_df[unplotted_repl, resp_idx]),
                          matrix(NA,
                                 ncol = 3,
                                 nrow = n_dose - length(resp_y_values)))
    unplt_repl_df <- data.frame(Cx_axis_values, y_unplt_repl, row.names = NULL)
    names(unplt_repl_df) <- c("Cx_axis_values", seq(length(unplotted_repl)))
    unplt_mdf <- reshape2::melt(unplt_repl_df,
                                id.vars = "Cx_axis_values",
                                variable.name = "Replicate")
    myplot <- myplot +
      geom_point(data = unplt_mdf, size = 2,
                 aes(x = Cx_axis_values, y = value, shape = Replicate)) +
      scale_shape(solid = TRUE) + scale_shape_manual(values = c(16, 17, 18))
    # get Scores
    score_matrix <- compute_mixpred_scores(mix_df,
                                           mix_idx,
                                           unplotted_repl,
                                           curve_data_list)
    score_record[iter_idx, ] <- c(mix_idx, array(t(score_matrix)))
    # to annotate plot, create strings of scores
    llh_vals <- sprintf("%.2f", score_matrix[1, ])
    mse_vals <- sprintf("%.2f", score_matrix[2, ])
    crps_vals <- sprintf("%.2f", score_matrix[3, ])
    summry <- paste(c("LLH:", llh_vals))
    summry2 <- paste(c("MSE:", mse_vals))
    summry3 <- paste(c("CRPS:", crps_vals))
    #for future use, may print scores on plots
    my_sums = c(summry, summry2, summry3) #nolint
    # annotation offset? set offset_y to 10% of max(curve_data_IA)
    # add chemical presence key
    active_chem <- (chem_conc_matr_reordered[1, ] > 0) + 1
    chem_cols <- c("white", "black")[active_chem]
    names(chem_cols) <- names(pure_ordering)
    # provide subtitle to plot with mixture description
    CAS_desc <- mix_guide$Description[which(mix_guide$CAS ==
                                              mix_df$CAS[mix_idx])]
    if (any(grep("Graded", CAS_desc, ignore.case = TRUE))) {
      # remove initial "Graded" text?
      CAS_desc <- sub("Graded [1-9] *", "", CAS_desc, ignore.case = TRUE)
      # remove double spaces
      CAS_desc <- sub(" {2,}", " ", CAS_desc, ignore.case = TRUE)
    }
    annotated_plot <- myplot +
      xlab("Concentration (uM)") +
      ylab("Response") + labs(title = mix_df$CAS[mix_idx],
                              subtitle = paste(strwrap(CAS_desc,
                                                       width = .7 *
                                                         getOption("width")),
                                               collapse = "\n"))
    # make some arbitrary plot and get the legend in
    fake_dats <- mtcars[1:18, 1:2]
    fake_dats$col <- names(pure_ordering)
    fake_plot <- ggplot2::ggplot(data = fake_dats, aes(mpg, cyl, color = col)) +
      geom_point(shape = 15) +
      scale_color_manual(name = "Active CAS", values = chem_cols)
    guide_color <- get_legend(fake_plot)
    # combine the response plot with the legend from the arbitrary plot using
    # command plot_grid
    print(plot_grid(annotated_plot + theme(legend.position = "right"),
                    guide_color,
                    rel_widths = c(.6, .2)))
  }
  return(score_record)
}

#' Compare individual and Mixture response
#'
#' Creates a series of plots where each plot shows the mixture response with one
#' of the individual component responses.  The mixture response concentration is
#' scaled to only correspond to the concentration of the compared chemical.  For
#' example, if a mixture has a concentraton of 5 and one component contribute 1,
#' when plotted against that component, the mixture response will be shifted to
#' a concentration of 1.
#'
#' This function is not optimized and has some hardcoded bounds.
#'
#' @param mix_idx integer specifying which mixture to plot
#' @param ymax positive real number for the plot y limit
#'
#' @return
#' @export
#'
#' @examples
plot_mix_vs_individuals <- function(mix_idx = 34, ymax = 120) {
  # explicit reference global variables
  y_i_T21 <- .GlobalEnv$y_i_T21
  resp_idx <- .GlobalEnv$resp_idx
  pure_unique_CAS <- .GlobalEnv$pure_unique_CAS
  # conc_mat is orderd by pure_unique
  chem_conc_matr <- get_conc_matrix(mix_idx)
  chem_conc_matr_reordered <- chem_conc_matr[, pure_ordering]
  active_idx <- which(chem_conc_matr_reordered[1, ] > 0)
  # plot mix compared to pure for each chemical
  par(mfrow = c(1, 2))
  if (length(active_idx) > 4) par(mfrow = c(3, 3))
  tot_df <- data.frame()
  for (i in active_idx) {
    chem_cas <- pure_unique_CAS[i]
    chem_conc <- chem_conc_matr_reordered[1:15, i]
    resp_y_values <- array(unlist(mix_df[mix_idx, resp_idx]))
    plot(chem_conc, resp_y_values,
         log = "x", type = "l", lwd = 3, col = 3,
         xlim = c(1e-12, 1e-4), ylim = c(0, ymax),
         main = paste(chem_cas, "vs Mix", sep = " "))
    matching_pure <- which(y_i_T21$CAS == chem_cas)
    y_i_unscaled <- y_i_T21[matching_pure, 3:ncol(y_i_T21)]
    matplot(x = t(Cx[matching_pure, ]),
            y = t(y_i_unscaled),
            type = "l", add = TRUE, col = 1, log = "x")
    n_repl_cas <- ncol(t(y_i_unscaled))
    df_part <- data.frame(chem_cas,
                          chem_conc,
                          resp_y_values,
                          Cx[matching_pure[1], ],
                          t(y_i_unscaled))
    names(df_part) <- c("CAS", "Mix Conc", "Mix Resp", paste("CAS Conc"),
                        paste("CAS Resp", 1:n_repl_cas))
    df_long <- reshape2::melt(df_part, id.vars = 1:4)
    tot_df <- rbind(tot_df, df_long)
  }
}

#' Visualize Feasible Clusters
#'
#' Creates a plot that visualizes the top clusters as connected line segments.
#'
#' @param slopes vector of real slope values
#' @param clust_centers_w_prob list object containing top clusters, returned
#'   from [cluster_centers()].
#'
#' @return original plot settings from command [par()]
#' @export
#'
#' @examples
visualize_clusters <- function(slopes, clust_centers_w_prob) {
  string2num <- function(xstr) as.numeric(unlist(strsplit(xstr, " ")))
  clust_names <- names(clust_centers_w_prob$assign)
  n_clust <- length(clust_names)
  oldpar <- par(mar = (c(5, 2, 4, 2) + .1))
  # hardcoded range, avoid?
  plot(
    x = slopes, y = rep(n_clust + 1, length(slopes)), xlim = c(min(slopes), 2),
    pch = 6, ylim = c(0, n_clust + 1), yaxt = "n", ylab = NA
  )
  for (clust_id in 1:length(clust_names)) {
    clust_assign <- string2num(clust_names[clust_id])
    for (cid in unique(clust_assign)) {
      active_x <- slopes[which(clust_assign == cid)]
      lines(x = c(range(active_x)), y = rep(n_clust - clust_id, 2), lwd = 2)
    }
  }
  par(oldpar)
}



#' Visualize Feasible Clusters
#'
#' Creates a plot that visualizes the top clusters as boxes surrounding the
#' point estimates.  Provides additional information about the ordering on the
#' left and the frequency of the particular clustering on the right.
#'
#' @param slopes a vector of the slope values
#' @param clust_centers_w_prob a named list with the results of the dirichlet
#'   process clustering
#' @param ul Upper limit for plotting the cluster values
#' @param ll lower limit for plotting the cluster values
#'
#' @return no return object specified other than original plot settings [par()]
#' @export
#'
#' @examples
visualize_clusters_blocks <- function(slopes,
                                      clust_centers_w_prob,
                                      ul = 2, ll = 0) {
  string2num <- function(xstr) as.numeric(unlist(strsplit(xstr, " ")))
  clust_names <- names(clust_centers_w_prob$assign)
  n_clust <- length(clust_names)
  oldpar <- par(mar = (c(5, 2, 4, 2) + .1))
  # hardcoded range, avoid?
  plot(
    x = slopes, y = rep(n_clust + 1, length(slopes)), xlim = c(ll, ul),
    pch = 1, ylim = c(0, n_clust - 1), yaxt = "n", ylab = NA,
    frame = FALSE, xlab = "Fitted Slope", main = "Most Likely Clusters"
  )
  axis(2, at = 0:n_clust, label = n_clust:0, las = 1, lwd = 0, lwd.ticks = 1)
  axis(4, at = 0:n_clust - 1,
       label = paste("(", clust_centers_w_prob$assign[n_clust:0 + 1], ")",
                     sep = ""),
       las = 1, lwd = 0, line = -1)
  for (clust_id in 1:length(clust_names)) {
    clust_assign <- string2num(clust_names[clust_id])
    for (cid in unique(clust_assign)) {
      points(
        x = slopes,
        y = rep(n_clust - clust_id,length(slopes)),
        xlim = c(min(slopes), 2),
        pch = ".", ylim = c(0, n_clust + 1), yaxt = "n", ylab = NA
      )
      active_x <- slopes[which(clust_assign == cid)]
      x_lims <- range(active_x) + c(-.01, +.01)
      y_lims <- rep(n_clust - clust_id, 2) + c(-.2, .2)
      rect(
        xleft = x_lims[1], xright = x_lims[2],
        ybottom = y_lims[1], ytop = y_lims[2]
      )
    }
  }
  par(oldpar)
}


# Violin Plot for Scores ####
plot_scores <- function(score_df, bootstrap_calc_list, method_levels = NA) {
  plot_df <- score_df
  names(plot_df) <- c("Mix Desc", rep(names(bootstrap_calc_list), 3))
  subdf <- plot_df # [subindx,]
  mdf <- rbind(cbind(reshape2::melt(subdf[, grep("CRPS", names(score_df))]),
                     "Score" = "CRPS"),
               cbind(reshape2::melt(subdf[, grep("MSE", names(score_df))]),
                     "Score" = "MSE"),
               cbind(reshape2::melt(subdf[, grep("LLH", names(score_df))]),
                     "Score" = "LLH"))
  names(mdf) <- c("Method", "Value", "Score")
  if (length(method_levels) > 1) {
    mdf$Method <- factor(mdf$Method, ordered = TRUE, levels = method_levels)
  }
  ggplot2::ggplot(mdf, aes(x = Method, y = Value)) +
    geom_boxplot(draw_quantiles = c(.25, .5, .75)) +
    facet_wrap(vars(Score), scales = "free_y") +
    scale_y_log10() +
    theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0)) +
    labs(title = "Score Summary")
}



plot_individual_response <- function(replicate_sets,
                                     re_par_summary,
                                     curve_fits,
                                     RE_curve_fits,
                                     Cx, y_i) {
  replicate_counts <- rep(1:18, unlist(lapply(replicate_sets, length)))
  # get MSE
  DRM_SSE <- 0
  REM_SSE <- 0
  #if saving plot: pdf(width=10, height = 8, onefile = FALSE)
  old_par <- par(mfrow = c(3, 3))
  for (ri in 1:length(replicate_sets)) {
    u_RE_vals <- re_par_summary$u_RE_params[which(replicate_counts == ri)]
    v_RE_vals <- re_par_summary$v_RE_params[which(replicate_counts == ri)]
    rep_idx <- replicate_sets[[ri]]
    plot(Cx[rep_idx, ], y_i[rep_idx, ],
         log = "x", pch = 1:length(rep_idx),
         xlab = "Concentration",
         ylab = "Response",
         main = chem_map[ri])
    lines(Cx[rep_idx[1], ],
          hill_function(curve_fits[ri, 1],
                        curve_fits[ri, 2],
                        curve_fits[ri, 3],
                        Cx[rep_idx[1], ]),
          col = 2)
    DRM_SSE <- DRM_SSE + sum((y_i[rep_idx, ] -
                                hill_function(curve_fits[ri, 1],
                                              curve_fits[ri, 2],
                                              curve_fits[ri, 3],
                                              Cx[rep_idx[1], ]))^2,
                             na.rm = TRUE)
    for (re_idx in 1:length(u_RE_vals)) {
      U_RE <- u_RE_vals[re_idx]
      V_RE <- v_RE_vals[re_idx]
      RE_fit_pred <- V_RE + hill_function(RE_curve_fits[ri, 1] + U_RE,
                                          RE_curve_fits[ri, 2],
                                          RE_curve_fits[ri, 3],
                                          Cx[rep_idx[re_idx], ])
      lines(Cx[rep_idx[re_idx], ], RE_fit_pred, col = 3)
      REM_SSE <- REM_SSE + sum((y_i[rep_idx[re_idx], ] - RE_fit_pred)^2,
                               na.rm = TRUE)
    }
  }
  #if saving plot: dev.off()
  total_points <- sum(!is.na(Cx))
  DRM_SSE / total_points
  REM_SSE / total_points
  # return plot settings
  par(old_par)
}



# Hill plots ####

# creates the Hill function example with reflected segments used in the
# manuscript
plot_reflection_illustration <- function(slope_in = .5) {
  a <- .9
  b <- .5
  slope_c <- slope_in
  # a: sill parameter
  # b:  EC50 or Disassociation constant
  # c: slope
  # d: minimum effect
  hilly_inverse_test <- function(y, c) hill_invs_factry(a, b, c)(y)
  y_test <- seq(a + .04, 2 * a - .01, by = sign(a) * .01)
  y_test_bey <- seq(2 * a - .01, 6 * a - .01, by = sign(a) * .01)
  y_test_neg <- seq(0, -4 * a - .01, by = -sign(a) * .01)
  x_test <- seq(0, 5, by = .01)
  lgnd_brks <- c(" GCA", "RGCA", "Extension", "Intermediate",
                 "Reflection", "Ext+Reflect", "alt")
  pt1 <- cbind(x_test, a / (1 + (b / x_test)^(slope_c)), 1, 1, lgnd_brks[2])
  pt1_ghost <- cbind(-x_test - 2 * b,
                     a / (1 + (b / x_test)^(slope_c)),
                     5, 1, lgnd_brks[4])
  pt2 <- cbind(sapply(y_test, function(x) hilly_inverse_test(x, slope_c)),
               y_test, 2, 1, lgnd_brks[5])
  pt3 <- cbind(sapply(y_test_neg, function(x) hilly_inverse_test(x, slope_c)),
               y_test_neg, 3, 1, lgnd_brks[3])
  pt4 <- cbind(sapply(y_test_bey, function(x) hilly_inverse_test(x, slope_c)),
               y_test_bey, 4, 1, lgnd_brks[6])
  x_test <- seq(-5, 5, by = .1)
  pt_gca <- cbind(x_test, a / (1 + (b / x_test)^(1)), 0, 2, lgnd_brks[1])
  pt_gca[which(pt_gca[, 2] == "Inf"), 2] <- NA
  hill_df <- data.frame(rbind(pt_gca, pt1, pt2, pt3, pt4, pt1_ghost))
  names(hill_df) <- c("x", "y", "subcurve", "line", "Segment")
  hill_df$x <- as.numeric(hill_df$x)
  hill_df$y <- as.numeric(hill_df$y)
  rect_region <- data.frame(xmin = c(0, -b),
                            xmax = c(5, 0),
                            ymin = c(0, -5),
                            ymax = c(a, 0),
                            z = c(1, 2))
  basic_col <- c("black", "green", "red", "gray", "cyan", "purple", "gray")
  gplot <- ggplot2::ggplot(data = hill_df,
                           aes(x = x, y = y, color = Segment)) +
    geom_vline(xintercept = 0, linetype = 1, color = "black") +
    geom_hline(yintercept = 0, linetype = 1, color = "black") +
    geom_rect(data = rect_region, inherit.aes = FALSE,
              aes(xmin = xmin,
                  xmax = xmax,
                  ymin = ymin,
                  ymax = ymax,
                  fill = as.factor(z)),
              alpha = .2,
              show.legend = FALSE) +
    geom_line(linewidth = 1,
              lineend = "round",
              linetype = as.numeric(hill_df$line)) +
    geom_hline(yintercept = a, linetype = 2, color = "gray") +
    annotate("text", x = 1, y = 1, label = "Sill", size = 3, col = "blue") +
    geom_vline(xintercept = -b, linetype = 2, color = "gray") +
    annotate("text",
             x = -.55,
             y = -1,
             label = "-EC50",
             angle = 90,
             size = 3,
             col = "blue") +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          legend.margin = margin()) +
    scale_colour_manual(values = basic_col, breaks = lgnd_brks) +
    xlab("Concentration (uM)") +
    coord_cartesian(xlim = c(-2.5, 2), ylim = c(-1.5, 3)) +
    guides(colour =
             guide_legend(nrow = 1, 
                          override.aes = 
                            list(linetype =c(3, 1, 1, 1, 1, 1)))) +
    ylab("Response") +
    labs(title = "Extending GCA with Reflection")
  x_test <- seq(0.01, 5, by = .01) # 1*10^(seq(-8, .5, by=.1))
  for (slope_c in c(seq(.1, .8, by = .05), seq(1, 4, by = .2))) {
    pt1 <- cbind(x_test, a / (1 + (b / x_test)^(slope_c)), 1, 1, "Alt")
    pt2 <- cbind(
      sapply(y_test, function(x) hilly_inverse_test(x, slope_c)),
      y_test, 2, 1, "Alt"
    )
    pt3 <- cbind(
      sapply(y_test_neg, function(x) hilly_inverse_test(x, slope_c)),
      y_test_neg, 3, 1, "Alt"
    )
    pt4 <- cbind(
      sapply(y_test_bey, function(x) hilly_inverse_test(x, slope_c)),
      y_test_bey, 4, 1, "Alt"
    )
    
    hill_df <- data.frame(rbind(pt_gca, pt1, pt2, pt3, pt4, pt1_ghost))
    names(hill_df) <- c("x", "y", "subcurve", "line", "Segment")
    hill_df$x <- as.numeric(hill_df$x)
    hill_df$y <- as.numeric(hill_df$y)
    # set last point of extension to NA to avoid complete line
    hill_df[which(hill_df$subcurve == "4")[1] - 1, ]$y <- NA
    gplot <- gplot +
      geom_line(
        data = hill_df,
        aes(x = x, y = y, color = Segment),
        linewidth = .25,
        linetype = 1
      )
  }
  # reorder layers
  n_layers <- length(gplot$layers)
  gplot$layers <- gplot$layers[c(1, 7:n_layers, 2:6)]
  gplot
}
# for saving, use pdf("RGCA_symmetry_full.pdf", width =7,height = 4)
# then plot_reflection_illustration(); dev.off()

# not used, but plots the inverse of the Hill function
plot_inverse_illustration <- function(plot_reciprocal = FALSE) {
  a <- -.9
  b <- .5
  slope_c <- 1.5
  # max_R allows for variable definition of partial agonist
  # rather than assuming partial means < 1 response
  # a: sill parameter
  # b:  EC50 or Disassociation constant
  # c: slope
  # d: minimum effect
  hilly_inverse_test <- function(y, c) hill_invs_factry(a, b, c)(y)
  #range   [a, 2a]
  y_test <- seq(a + .04, 2 * a - .01, by = sign(a) * .01)
  #range  [2a, Inf (~6a)]
  y_test_bey <- seq(2 * a - .01, 6 * a - .01, by = sign(a) * .01)
  y_test_neg <- seq(0, -4 * a - .01, by = -sign(a) * .01) # [0, -4a]
  x_test <- seq(0, 5, by = .01)
  lgnd_brks <- c(" GCA",
                 "RGCA",
                 "Extension",
                 "Intermediate",
                 "Reflection",
                 "Ext+Reflect")
  pt1 <- cbind(a / (1 + (b / (x_test))^(slope_c)), x_test, 1, 1, lgnd_brks[2])
  pt2 <- cbind(y_test,
               sapply(y_test, function(x) hilly_inverse_test(x, slope_c)),
               2, 1, lgnd_brks[5])
  pt3 <- cbind(y_test_neg,
               sapply(y_test_neg, function(x) hilly_inverse_test(x, slope_c)),
               3, 1, lgnd_brks[3])
  pt4 <- cbind(y_test_bey,
               sapply(y_test_bey, function(x) hilly_inverse_test(x, slope_c)),
               4, 1, lgnd_brks[6])
  x_test <- seq(-5, 5, by = .01)
  pt_gca <- cbind(a / (1 + (b / x_test)^(1)), x_test, 0, 2, lgnd_brks[1])
  hill_df <- data.frame(rbind(pt_gca, pt1, pt2, pt3, pt4))
  names(hill_df) <- c("x", "y", "subcurve", "line", "Segment")
  hill_df$x <- as.numeric(hill_df$x)
  hill_df$y <- as.numeric(hill_df$y)
  if (plot_reciprocal) {
    hill_df$y <- 1 / as.numeric(hill_df$y)
  }
  hill_df$y[which(hill_df$y == "Inf")] <- NA
  basic_col <- c("#E69F00",
                 "#56B4E9",
                 "#009E73",
                 "#F0E442",
                 "#0072B2",
                 "#D55E00",
                 "#CC79A7")
  basic_col <- c("black", "green", "red", "gray", "cyan", "purple")
  gplot <- ggplot2::ggplot(data = hill_df, aes(x = x, y = y, color = Segment)) +
    geom_line(linewidth = 1,
              lineend = "round",
              linetype = as.numeric(hill_df$line)) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          legend.margin = margin()) +
    scale_colour_manual(values = basic_col, breaks = lgnd_brks) +
    xlab("Concentration (uM)") +
    guides(colour = guide_legend(nrow = 1,
                                 override.aes =
                                   list(linetype = c(3, 1, 1, 1, 1)))) +
    ylab("Response") +
    labs(title = "Extending GCA with Reflection, Inverse")
  if (plot_reciprocal) {
    gplot <- gplot + coord_cartesian(xlim = c(-1, 5), ylim = c(-5, 5))
  } else {
    gplot <- gplot + coord_cartesian(xlim = c(-3, 4), ylim = c(-3, 4))
  }
  print(gplot)
}

# No Effect Inclusion Test #####
# Include/exclude no effects testing:
# for a mixture with chemicals that have 0 effect, does their
# inclusion affect the mixture response?
compare_exclude_include <- function() {
  plot_compare <- function(idx_set_list, conc_factor = 1, plt_title = NA) {
    reptit <- sapply(1:length(idx_set_list),
                     FUN = function(x) length(idx_set_list[[x]]))
    EC_df <- cbind(mix_df[unlist(idx_set_list), ],
                   "group" = rep(names(idx_set_list), reptit),
                   "multp" = rep(conc_factor, reptit))
    # calibrate actual dosage using the guide
    plot_df <- cbind(EC_df[, c(2, 4, 6)], EC_df[, c(6 + 15, 51, 52)])
    for (idxoff in 1:14) {
      sub_df <- cbind(EC_df[, c(2, 4, 6 + idxoff)],
                      EC_df[, c(6 + idxoff + 15, 51, 52)])
      names(sub_df) <- names(plot_df)
      plot_df <- rbind(plot_df, sub_df)
    }
    plot_df$ReplicateSet <- as.factor(plot_df$ReplicateSet)
    plot_df$CONC0 <- plot_df$CONC0 * plot_df$multp
    p <- ggplot2::ggplot(plot_df, aes(x = CONC0,
                                      y = RESP0,
                                      group = group,
                                      color = group,
                                      shape = ReplicateSet)) +
      geom_point(size = 3, stroke = 1.5) +
      scale_x_continuous(trans = "log10") +
      xlab("Concentration (uM)") +
      ylab("Response") +
      labs(title = plt_title)
    return(p)
  }
  EC_AR_guide <- as.numeric(mix_guide[which(mix_guide$CAS == "NOCAS_49497"),
                                      4:21])
  EC_all_guide <- as.numeric(mix_guide[which(mix_guide$CAS == "NOCAS_49538"),
                                       4:21])
  EC_ER_guide <- as.numeric(mix_guide[which(mix_guide$CAS == "NOCAS_49517"),
                                      4:21])
  crat_1 <- EC_AR_guide / EC_all_guide
  mean(crat_1[crat_1 > 0])
  crat_2 <- EC_ER_guide / EC_all_guide
  mean(crat_2[crat_2 > 0])
  idx_ec_ar <- which(mix_df$CAS == "NOCAS_49497") # EquiConc AR
  idx_ec_all <- which(mix_df$CAS == "NOCAS_49538") # EquiConc all
  idx_ec_er <- which(mix_df$CAS == "NOCAS_49517") # EquiConc ER
  idx_ep_ar2 <- which(mix_df$CAS == "NOCAS_49561") # EquiConc AR alt
  idx_set_list <- list("EC_all" = idx_ec_all,
                       "EC_AR" = idx_ec_ar,
                       "EC_ER" = idx_ec_er,
                       "EC_AR2" = idx_ep_ar2)
  p1 <- plot_compare(idx_set_list,
                     conc_factor = c(1, 1, 3, 1),
                     plt_title = "Equiconcentration Comparison")
  EP_all_guide <- as.numeric(mix_guide[which(mix_guide$CAS == "NOCAS_49537"),
                                       4:21])
  EP_AR_guide <- as.numeric(mix_guide[which(mix_guide$CAS == "NOCAS_49496"),
                                      4:21])
  EP_ER_guide <- as.numeric(mix_guide[which(mix_guide$CAS == "NOCAS_49516"),
                                      4:21])
  crat_1 <- EP_AR_guide / EP_all_guide
  mean(crat_1[crat_1 > 0])
  crat_2 <- EP_ER_guide / EP_all_guide
  mean(crat_2[crat_2 > 0])
  idx_ep_all <- which(mix_df$CAS == "NOCAS_49537") # EquiPot all
  idx_ep_ar <- which(mix_df$CAS == "NOCAS_49496") # EquiPot AR
  idx_ep_er <- which(mix_df$CAS == "NOCAS_49516") # EquiPot ER
  EP_set_list <- list("EP_all" = idx_ep_all,
                      "EP_AR" = idx_ep_ar,
                      "EP_ER" = idx_ep_er)
  p2 <- plot_compare(EP_set_list,
                     conc_factor = c(1, 4, 2),
                     plt_title = "Equipotent Comparison")
  plot_grid(p1, p2)
  # EC_all and EC_AR may be mixed up!
}

# Sim Study ####
sim_study_boxplot <- function(record_MSE, data_col_names, record_CRPS = NA) {
  # Plot results as boxplot
  if (!is.na(record_CRPS)) {
    record <- data.frame(rbind(cbind(record_MSE, "Score" = "MSE"),
                               cbind(record_CRPS, "Score" = "CRPS")))
    names(record) <- data_col_names
    record_df <- reshape2::melt(record, id.vars = c("Score", "Mix_ID"))
    names(record_df) <- c("Score_Method", "Mix_ID", "Pred_Method", "Score")
    record_df$Score <- as.numeric(record_df$Score)
    my_gplot <- ggplot2::ggplot(record_df,
                                aes(x = as.factor(Mix_ID), y = Score)) +
      geom_boxplot(aes(colour = Pred_Method)) +
      facet_grid(Score_Method ~ Mix_ID, scales = "free") +
      theme_bw() +
      xlab("True Mixture Model") +
      theme(legend.position = "bottom") +
      guides(colour = guide_legend(title = "Fitting Method"))
    return(my_gplot)
  }
  record <- data.frame(record_MSE)
  names(record) <- data_col_names
  record_df <- reshape2::melt(record, id.vars = c("Num Chems", "Mix ID"))
  names(record_df) <- c("Num Chems", "Mix ID", "Pred_Method", "MSE")
  record_df$MSE <- as.numeric(record_df$MSE)
  record_df$`Num Chems` <- as.numeric(record_df$`Num Chems`)
  record_df <- record_df[record_df$MSE < 20, ]
  ggplot2::ggplot(record_df, aes(x = as.factor(`Mix ID`), y = MSE)) +
    geom_boxplot(aes(colour = Pred_Method)) +
    facet_grid(`Num Chems` ~ ., scales = "free") +
    theme_bw() +
    xlab("True Mixture Model") +
    theme(legend.position = "bottom") +
    coord_cartesian(ylim = c(0, 4)) +
    guides(colour = guide_legend(title = "Fitting Method"))
}

# Isobole ####
create_RGCA_isobole_plot <- function() {
  isobole_df <- data.frame()
  for (a_val in c(-2, 1, 3, 5)) {
    for (slope_val in c(.5, 1, 1.5)) {
      # create a parameter input matrix
      param_matrix <- as.matrix(cbind("a" = c(a_val, 5),
                                      "b" = c(.7, 1),
                                      "c" = c(slope_val, 1),
                                      "max_R" = 5,
                                      "d" = 0))
      # specify both chems in cluster 1 of 1
      cluster_assign_vec <- c(1, 1)
      # create a calculator to predict response given concentration
      mix_pred <- mix_function_generator(param_matrix, cluster_assign_vec)
      # create matrix of concentration values
      c_val <- seq(.05, 5, by = .01)
      n_pix <- length(c_val)
      c_mat <- expand.grid(c_val, c_val)
      # apply calculator to concentrations
      pred_effect <- apply(c_mat, MARGIN = 1, FUN = function(x) mix_pred(x))
      mix_df <- data.frame(cbind(c_mat, pred_effect))
      names(mix_df) <- c("C1", "C2", "effect")
      isobole_df <- rbind(isobole_df,
                          cbind(mix_df, "a" = a_val, "slope" = slope_val))
    }
  }
  # 2d contour plot
  #save with pdf(file='Isoboles_RGCA_x.pdf',width=8, height = 6)
  ggplot2::ggplot(isobole_df, aes(x = C1, y = C2, z = effect)) +
    geom_contour(bins = 30) +
    facet_grid(cols = vars(a), rows = vars(slope))
  #close save file dev.off()
  # 3d plot with 2d contour plot on floor
  plotly::plot_ly(x = (c_val),
                  y = (c_val),
                  z = matrix(mix_df$effect, nrow = n_pix)) %>%
    plotly::add_surface(contours = list(z = list(show = TRUE,
                                                 usecolormap = TRUE,
                                                 highlightcolor = "#ff0000",
                                                 project = list(z = TRUE)),
                                        zlim = c(0, 3)))
  # can:  %>% layout(scene = list(camera=list(eye = list(x=2, y=.5, z=.7))))
  # can:  %>% add_surface(z = matrix(pred_regular, nrow = n_pix), opacity = .4)
  # can:  %>% add_surface(z = matrix(pred_agonist, nrow = n_pix), opacity = .4)
  # alternate: plot_ly(z=matrix(mix_df$effect, nrow = n_pix), type = "contour")
}

# Uniqueness and Stability ####
plot_RGCA_soln_intersect <- function() {
  solutions_df <- data.frame()
  # test: # -0.575: one solution .45, -.35,-.01, , -.75, -.35, -.1
  for (a_val in c(-2, -.4, -.1)) {
    for (slope_val in c(1.5, 1, 15)) {
      for (c1 in 3) {
        for (c2 in 1) {
          # if using a third chem: c3 <- 1
          param_matrix <- as.matrix(cbind("a" = c(a_val, 1),
                                          "b" = c(1, .6),
                                          "c" = c(slope_val, 1),
                                          "max_R" = 1))
          hill_inverse_list <- apply(param_matrix,
                                     MARGIN = 1,
                                     function(x) {
                                       do.call(hill_invs_factry, as.list(x))
                                     })
          myfun <- function(r) {
            c1 / hill_inverse_list[[1]](r) +
              c2 / hill_inverse_list[[2]](r)
            # if additional chem: c3/hill_inverse_list[[3]](r)
          }
          myseq <- seq(-1, 1, by = .001)
          curr_df <- data.frame(rbind(
            cbind(myseq, c1 / sapply(myseq, hill_inverse_list[[1]]), "C1", 1),
            cbind(myseq, c2 / sapply(myseq, hill_inverse_list[[2]]), "C2", 2),
            #add cbind(myseq, c3/sapply(..)) if there is third chem
            cbind(myseq, sapply(myseq, myfun), "C3", 4)
          ))
          names(curr_df) <- c("R", "x", "curveID", "color")
          solutions_df <- rbind(solutions_df, cbind(curr_df,
                                                    "sill" = a_val,
                                                    "slope" = slope_val,
                                                    "c1" = c1,
                                                    "c2" = c2))
        }
      }
    }
  }
  solutions_df$R <- as.numeric(solutions_df$R)
  solutions_df$x <- as.numeric(solutions_df$x)
  solutions_df$c2 <- as.factor(solutions_df$c2)
  solutions_df$sill <- paste("Sill:", solutions_df$sill)
  solutions_df$slope <- paste("Slope:", solutions_df$slope)
  ggplot2::ggplot(solutions_df, aes(x = R, y = x, color = color)) +
    geom_line(linewidth = 1) +
    ylim(-20, 20) +
    geom_hline(yintercept = 1, linetype = 2) +
    facet_grid(sill ~ slope) +
    scale_colour_viridis_d("Chemical",
                           labels = c("Chem1",
                                      "Chem2",
                                      "Chem3",
                                      "Mix"), direction = -1) +
    theme(legend.position = "bottom")
  #manual save with ggsave("RGCA_limit_-2_1p5.pdf", device = "pdf", 8x9 in)
}

# MCMC chain analysis
plot_MCMC_chains <- function(re_chains) {
  burnin <- 5000
  if (re_iter < burnin) burnin <- 100
  thin_idx <- seq(burnin, re_iter - 100, by = 20)
  # paroi:  parameter of interest
  paroi <- re_chains$sill_mideffect_record[thin_idx, 17 + 0 * n_chems]
  points(paroi, log = "y") # , ylim = c(1e-8, 1000))
  paroi <- re_chains$sill_mideffect_record[thin_idx, 17 + 0 * n_chems]
  paroi_sl <- re_chains$slope_record[thin_idx, 17]
  paroi_ec <- re_chains$sill_mideffect_record[thin_idx, 17 + 1 * n_chems]
  plot(paroi)
  plot(paroi_ec, log = "y")
  plot(paroi_sl, log = "y")
  plot(re_chains$slope_record[thin_idx, 1] /
         ((re_chains$sill_mideffect_record[thin_idx, 1 + n_chems]) /
            (re_chains$sill_mideffect_record[thin_idx, 1])), log = "y")
  plot(log(re_chains$tau[thin_idx]))
  # compare:  6 vs 1 vs 17
  #if saving, pdf("MCMC_traces.pdf", width = 10, height = 8)
  par(mfrow = c(3, 4))
  for (idx in c(1, 6, 17)) {
    paroi <- re_chains$sill_mideffect_record[thin_idx, idx + 0 * n_chems]
    paroi_sl <- re_chains$slope_record[thin_idx, idx]
    paroi_ec <- re_chains$sill_mideffect_record[thin_idx, idx + 1 * n_chems]
    plot(re_chains$sigma[thin_idx, idx])
    # if(idx == 17){plot(paroi, ylim = c(-.5, 7),
    # main = paste("Sill Trace plot, chemical ", idx), ylab = NA)
    # }else
    plot(paroi, main = paste("Sill Trace plot, chemical ", idx), ylab = NA)
    plot(paroi_ec, log = "y",
         main = paste("EC50 Trace plot (log y), chemical ", idx),
         ylab = NA)
    plot(paroi_sl, main = paste("Slope Trace plot, chemical ", idx), ylab = NA)
  }
  #close save file  dev.off()
  plot(re_chains$sill_mideffect_record[thin_idx, idx + 1 * n_chems], log = "y")
  plot(re_chains$slope_record[thin_idx, idx], log = "y")
  thin_idx <- seq(10000, re_iter - 100, by = 40)
  plot(hill_function(
    re_chains$sill_mideffect_record[thin_idx, idx + 0 * n_chems],
    re_chains$sill_mideffect_record[thin_idx, idx + 1 * n_chems],
    re_chains$slope_record[thin_idx, idx],
    1e-10
  ))
}
