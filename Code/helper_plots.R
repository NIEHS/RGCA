## Plotting functions for DP+RGCA method

library(ggplot2)
library(cowplot)
library(reshape2)
library(scoringRules)

#' Plot individual dose response curves
#'
#' Plots a few hard-coded samples of the data for a comparison used in the
#' manuscript.  Shows why a random effect model is needed for the Tox21 data.
#' Requires the data cleaning script to be previously sourced,
#' source("Code/tox21_prep_data.R").
#'
#' @param df a data frame with the prepared tox21 data from source("Code/tox21_prep_data.R")
#'
#' @return a ggplot with four examples of the tox21 data
#' @export
#'
#' @examples
plot_sample_indv = function(df){
  chems_to_plot = c("Hydroxyflutamide", "Progesterone" ,"Fluoxymestrone","Ethylenediamine")
  # get the names of all chemicals in the data
  sample_names = unique(pure_df$Sample.Name)
  #extract x-y values and a map of replicate to idx
  x_vals = c()
  y_vals = c()
  
  for(chem_nm in chems_to_plot){
    cm_idx = which(df$Sample.Name == chem_nm)
    # add CAS to chem label for clarity
    chem_CAS = names(chem_map_plain)[which(chem_map_plain == chem_nm)]
    chem_label = paste(chem_nm, " (",chem_CAS,")", sep = "")
    # create long-format data frame with concentrations (x-values)
    new_x = cbind(melt(t(df[cm_idx, conc_idx])),"chem"= chem_label)
    # substitute replicate index with replicate number
    new_x$Var2 = sapply(new_x$Var2,FUN =function(x) which(x == cm_idx))
    # create long-format data with responses (y-value)
    new_y = cbind(melt(t(df[cm_idx, resp_idx])))
    x_vals = rbind(x_vals, new_x)
    y_vals = rbind(y_vals, new_y)
  }
  # combine into a data frame and subset to drop duplicates
  gg_df = cbind(x_vals, y_vals)[,c(4,2,3,7)]
  names(gg_df) = c("chem","idx","Concentration", "Response")
  ggplot(data = gg_df, aes(Concentration, Response, colour = as.factor(idx) ))+
    geom_point()+geom_line()+scale_x_continuous(trans = "log10")+
    theme(legend.position = "top")+
    facet_wrap(vars(chem), scales = "free")+
    guides(colour = guide_legend(nrow = 1))+
    labs(colour='Tox21 AR-luc Replicates')
  
  
}
#pdf("Output/replicate_ex.pdf", width =7, height = 4)
#plot_sample_indv(df)
#dev.off()



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
#' 
make_ggplot = function(Cx_axis_values,resp_y_values, curve_data, curve_data_GCA = NULL, curve_data_IA = NULL){
  curve_data_quantiles = t(apply(curve_data, MARGIN = 2,
                                 FUN = function(x) quantile(x, probs = c(0.05, .5, .95), na.rm=T)))
  gg_df = data.frame(Cx_axis_values,resp_y_values, curve_data_quantiles)
  names(gg_df) = c("Cx", "resp_y", "low", "med", "up")
  hill_pred_plot = ggplot(gg_df, aes(Cx, resp_y))+
    geom_point(shape = 21, colour = "black", fill = "white", size = 3, stroke=1.5)+
    scale_x_continuous(trans='log10')+
    geom_line(aes(Cx, med, colour = "RE+DP"))+
    geom_ribbon(aes(ymin = low, ymax = up, fill = "RE+DP"), alpha=0.1, #fill = "green", 
                color = "black", linetype = "dotted")
  if(!is.null(curve_data_GCA)){
    curve_GCA_data_quantiles = t(apply(curve_data_GCA, MARGIN = 2,
                                       FUN = function(x) quantile(x, probs = c(0.05, .5, .95), na.rm=T)))
    GCA_df = data.frame(Cx_axis_values, resp_y_values, curve_GCA_data_quantiles)
    names(GCA_df) = c("Cx", "resp_y", "GCA_low", "GCA_med", "GCA_up")
    hill_pred_plot = hill_pred_plot+  
      geom_line(data=GCA_df, aes(Cx, GCA_med,colour = "GCA"))+
      geom_ribbon(data=GCA_df, aes(ymin = GCA_low, ymax = GCA_up, fill= "GCA"), 
                  alpha=0.1, #fill = "red", 
                  color = "black", linetype = "dotted")+
      scale_colour_manual(name = "Legend",values = c("RE+DP" = "black", "GCA" ="red"))+
      scale_fill_manual(name = "Legend", values = c("RE+DP" = "green", "GCA" ="red"))
  }
  if(!is.null(curve_data_IA)){
    curve_IA_data_quantiles = t(apply(curve_data_IA, MARGIN = 2,
                                      FUN = function(x) quantile(x, probs = c(0.05, .5, .95), na.rm=T)))
    IA_df = data.frame(Cx_axis_values, resp_y_values, curve_IA_data_quantiles)
    names(IA_df) = c("Cx", "resp_y", "IA_low", "IA_med", "IA_up")
    hill_pred_plot = hill_pred_plot+  
      geom_line(data=IA_df, aes(Cx, IA_med,colour = "IA"))+
      geom_ribbon(data=IA_df, aes(ymin = IA_low, ymax = IA_up, fill= "IA"), 
                  alpha=0.1, 
                  color = "black", linetype = "dotted")+
      scale_colour_manual(name = "Legend",values = c("RE+DP" = "black", "GCA" ="red", "IA" ="blue"))+
      scale_fill_manual(name = "Legend", values = c("RE+DP" = "green", "GCA" ="red", "IA" = "blue"))
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
make_ggplot_many = function(Cx_axis_values,resp_y_values, curve_data_list){
  # plt_id specified from the name of the curve calculator
  plt_id = names(curve_data_list)
  num_curves = length(curve_data_list)
  if(num_curves==0) return(NA)
  if(length(plt_id)!= num_curves) print("curves not named")
  
  quantile_df_colnames = c("Cx", "resp_y", "low", "med", "up", "ID")
  
  curve_data_quantiles = t(apply(curve_data_list[[1]], MARGIN = 2,
                                 FUN = function(x) quantile(x, probs = c(0.05, .5, .95), na.rm=T)))
  gg_df = data.frame(Cx_axis_values,resp_y_values, curve_data_quantiles, plt_id[1])
  names(gg_df) = quantile_df_colnames
  hill_pred_plot = ggplot(gg_df, aes(Cx, resp_y))+
    geom_point(shape = 21, colour = "black", fill = "white", size = 3, stroke=1.5)+
    scale_x_continuous(trans='log10')+
    geom_line(aes(Cx, med, colour = ID))+
    geom_ribbon(aes(ymin = low, ymax = up, fill = plt_id[1]), alpha=0.1, #fill = "green", 
                color = "black", linetype = "dotted")
  
  for(curve_idx in 2:num_curves){
    
    curve_data_quantiles = t(apply(curve_data_list[[curve_idx]], MARGIN = 2,
                                   FUN = function(x) quantile(x, probs = c(0.05, .5, .95), na.rm=T)))
    gg_df = data.frame(Cx_axis_values, resp_y_values, curve_data_quantiles,plt_id[curve_idx] )
    names(gg_df) = quantile_df_colnames
    hill_pred_plot = hill_pred_plot+  
      geom_line(data=gg_df, aes(Cx, med,colour = ID))+
      geom_ribbon(data=gg_df, aes(ymin = low, ymax = up, fill= ID), 
                  alpha=0.1, #fill = "red", 
                  color = "black", linetype = "dotted")
  }
  # coloring:  assume first plot is our method, use black line with green curve
  line_colors <- c("black", "red", "blue","orange", "magenta",  "yellow")
  fill_colors <- c("green", "red", "blue","orange", "magenta",  "yellow")
  
  
  
  legend_color_values =  1:num_curves
  fill_vals = 1:num_curves
  fill_vals[1]=3
  if(num_curves>2){
    legend_color_values = 1:(num_curves+1)
    legend_color_values=legend_color_values[-3]
    fill_vals = 1:(num_curves+1)
    fill_vals[1]= 3
    fill_vals=fill_vals[-3]
  }
  if(num_curves<length(line_colors)){
    legend_color_values =  line_colors[1:num_curves]
    fill_vals = fill_colors[1:num_curves]
  }
  
  
  names(legend_color_values) = plt_id
  names(fill_vals) = plt_id
  hill_pred_plot=hill_pred_plot+
    scale_colour_manual(name = "Legend",values =legend_color_values)+
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
#' @param test_idx
#'
#' @return
#' @export
#'
#' @examples
plot_dummy_mixture = function(Cx, y_i, tot_par_list, replicate_sets, 
                              bootstrapped_calc_list, 
                              test_idx = 1){
  #test along original chemical concentrations (equal dose)
  n_dose = ncol(Cx)
  dummy_conc_mat = matrix(0, nrow = n_dose, ncol = n_chems)
  dummy_conc_mat[, test_idx] = Cx[test_idx,]
  # curve list = c(RE_pred, GCA_pred, IA_pred)
  # curve_list = predict_mix_response(sampled_mix_funs,sampled_mix_funs_GCA,  
  #                                   sampled_mix_funs_IA, n_dose,
  #                                   dummy_conc_mat, default_entry=NA)
  curve_list = predict_mix_response_many(n_dose = n_dose, chem_conc_matr = dummy_conc_mat,
                                         bootstrap_calc_list = bootstrap_calc_list, 
                                         default_entry = NA)
  fit_resp = hill_function(tot_par_list$sill_params[test_idx], 
                           tot_par_list$ec50_params[test_idx],
                           tot_par_list$centers[test_idx], 
                           Cx[test_idx,])
  specific_idx = replicate_sets[[test_idx]][1]
  obs_resp = y_i[specific_idx,]
  #make_ggplot(Cx[specific_idx,],obs_resp,curve_list[[1]],curve_list[[2]], curve_list[[3]])
  make_ggplot_many(Cx[specific_idx,],obs_resp, curve_list)
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
plot_mixture_response = function(plot_set, mix_df, mix_conc_df, mix_guide, 
                                 bootstrap_calc_list ){
  n_calcs = length(bootstrap_calc_list)
  # preallocate matrix for score record
  score_record = matrix(0, nrow = sum(mix_df[plot_set,]$ReplicateSet==1), ncol = 1+n_calcs*3)
  for(iter_idx in 1:length(plot_set)){
    mix_idx = plot_set[iter_idx]
    chem_conc_matr = get_conc_matrix(mix_idx)
    n_dose = nrow(chem_conc_matr)-10 # first 15 are real, second 15 are extrapolation.  Here, just extrapolating lesss
    chem_conc_matr_reordered = chem_conc_matr[,pure_ordering]
    # curve_data_list = predict_mix_response(sampled_mix_funs,sampled_mix_funs_GCA,  
    #                                        sampled_mix_funs_IA, n_dose,
    #                                        chem_conc_matr_reordered)
    curve_data_list = predict_mix_response_many(n_dose,chem_conc_matr_reordered,
                                                bootstrap_calc_list)
    #Cx_axis_values = array(unlist(mix_df[mix_idx,conc_idx_T21_matrix]))
    Cx_axis_values = array(unlist(mix_conc_df[mix_idx,]))[1:n_dose]
    mix_replicates = which(mix_df$CAS ==mix_df$CAS[mix_idx] )
    resp_y_values = array(unlist(mix_df[mix_idx,resp_idx]))
    # allow extrapolation; match dim of Cx by adding null values to y observed
    y_response_ext = c(resp_y_values, rep(NA, n_dose-length(resp_y_values)))
    myplot = make_ggplot_many(Cx_axis_values, 
                              y_response_ext,
                              curve_data_list)
    
    # add replicates to single plot
    other_repl_idx = which(mix_df$CAS == mix_df$CAS[mix_idx])
    unplotted_repl  = setdiff(other_repl_idx, mix_idx)
    y_unplt_repl = rbind(t(mix_df[unplotted_repl, resp_idx]), matrix(NA,ncol=2, nrow=n_dose - length(resp_y_values)))
    unplt_repl_df = data.frame(Cx_axis_values, y_unplt_repl, row.names = NULL)
    names(unplt_repl_df) = c("Cx_axis_values", 1:length(unplotted_repl)+1)
    unplt_mdf = melt(unplt_repl_df, id.vars = "Cx_axis_values", variable.name = "Replicate")
    myplot = myplot+geom_point(data = unplt_mdf,
                               aes(x=Cx_axis_values, y=value,shape = Replicate))
    
    # add plot showing chemicals included
    
    
    # get Scores
    score_matrix = compute_mixpred_scores(mix_df, mix_idx, unplotted_repl, curve_data_list )
    score_record[iter_idx,] =  c(mix_idx,array(t(score_matrix)))
    # to annotate plot, create strings of scores
    llh_vals = sprintf("%.2f", score_matrix[1,])
    mse_vals = sprintf("%.2f", score_matrix[2,])
    crps_vals = sprintf("%.2f",score_matrix[3,])
    summry = paste(c("LLH:",llh_vals))
    summry2 = paste(c("MSE:",mse_vals))
    summry3 = paste(c("CRPS:", crps_vals))
    # annotation offset: 10% of max(curve_data_IA)
    #offset_y = max(unlist(curve_data_list))*.05
    
    # add chemical presence key
    active_chem = (chem_conc_matr_reordered[1,]>0)+1
    chem_cols = c("white", "black")[active_chem]
    names(chem_cols) = names(pure_ordering)
    
    # provide subtitle to plot with mixture description
    CAS_desc = mix_guide$Description[ which(mix_guide$CAS == mix_df$CAS[mix_idx])]
    if(any(grep("Graded", CAS_desc, ignore.case = TRUE))){
      # TODO: remove initial "Graded" text?
      CAS_desc = sub("Graded [1-9] *", "", CAS_desc, ignore.case = TRUE)
      # remove double spaces
      CAS_desc = sub(" {2,}", " ", CAS_desc, ignore.case = TRUE)
    }
    x_loc_ann = 5e-5 #5e-8
    annotated_plot = myplot+
      xlab("Concentration (uM)")+
      ylab("Response")+labs(title = mix_df$CAS[mix_idx], 
                            subtitle = paste(strwrap(CAS_desc, width = 1*getOption("width")), collapse = "\n"))
    # annotate("text", x=x_loc_ann, y=75+offset_y, label= summry3)+
    # annotate("text", x=x_loc_ann, y=75, label= summry2)+
    # annotate("text", x=x_loc_ann, y=75-offset_y, label= "(RE, GCA, IA)")+
    # #print(annotated_plot+scale_color_manual(name='Active Chemicals',values=chem_cols))
    
    # make some arbitrary plot and get the legend in
    #scm_plot
    fake_dats = mtcars[1:18,1:2]
    fake_dats$col =names(pure_ordering)
    fake_plot = ggplot(data = fake_dats,aes(mpg, cyl, color = col))+geom_point(shape = 15)+
      scale_color_manual(name='Active CAS', values=chem_cols)
    guide_color <- get_legend(fake_plot)
    
    # combine the response plot with the legend from the arbitrary plot using
    # plot_grid
    print(plot_grid(annotated_plot+theme(legend.position = "right"),
                    guide_color,rel_widths = c(.6, .2)))
  }
  return(score_record)
}


# compare ind vs mix,  mix # 34 as example for antagonism:
# pure #6: Oxymetholone "434-07-1"
# pure #9:  Hydroxyflutamide "52806-53-8"


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
plot_mix_vs_individuals = function(mix_idx = 34, ymax = 120){
  # conc_mat is orderd by pure_unique
  chem_conc_matr = get_conc_matrix(mix_idx)
  n_dose = nrow(chem_conc_matr)-10 # first 15 are real, second 15 are extrapolation.  Here, just extrapolating lesss
  chem_conc_matr_reordered = chem_conc_matr[,pure_ordering]
  active_idx = which(chem_conc_matr_reordered[1,]>0)
  #plot mix compared to pure for each chemical
  par(mfrow=c(1,2))
  if(length(active_idx)>4) par(mfrow=c(3,3))
  tot_df =data.frame()
  for(i in active_idx){
    chem_cas = pure_unique_CAS[i]
    chem_conc = chem_conc_matr_reordered[1:15,i]
    resp_y_values = array(unlist(mix_df[mix_idx,resp_idx]))
    plot(chem_conc, resp_y_values, log="x", type = "l",  lwd = 3, col=3,
         xlim = c(1e-12, 1e-4), ylim = c(0, ymax), 
         main = paste(chem_cas, "vs Mix", sep=" "))
    matching_pure = which(y_i_T21$CAS == chem_cas )
    y_i_unscaled  = y_i_T21[matching_pure,3:ncol(y_i_T21)]
    matplot(x=t(Cx[matching_pure,]), y=t(y_i_unscaled), 
            type = "l", add=T, col=1, log="x")
    n_repl_cas = ncol(t(y_i_unscaled)) 
    
    df_part = data.frame(chem_cas, chem_conc, resp_y_values, Cx[matching_pure[1],],t(y_i_unscaled))
    names(df_part)= c("CAS", "Mix Conc", "Mix Resp", 
                      paste("CAS Conc"),
                      paste("CAS Resp", 1:n_repl_cas))
    df_long = reshape2::melt(df_part,id.vars=1:4 )
    tot_df = rbind(tot_df, df_long)
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
visualize_clusters = function(slopes, clust_centers_w_prob){
  string2num = function(xstr) as.numeric(unlist(strsplit(xstr, " "))) 
  clust_names = names(clust_centers_w_prob$assign)
  n_clust = length(clust_names)
  oldpar = par(mar = (c(5,2,4,2)+.1))
  # hardcoded range, avoid?
  plot(x=slopes, y=rep(n_clust+1, length(slopes)), xlim =c(min(slopes), 2),
       pch = 6, ylim = c(0,n_clust+1), yaxt = "n", ylab = NA)
  for(clust_id in 1: length(clust_names)){
    clust_assign = string2num(clust_names[clust_id])
    for( cid in unique(clust_assign)){
      active_x = slopes[which(clust_assign==cid)]
      lines(x = c(range(active_x)), y=rep(n_clust-clust_id,2), lwd = 2)
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
#'
#' @return no return object specified other than original plot settings [par()]
#' @export
#'
#' @examples
visualize_clusters_blocks = function(slopes, clust_centers_w_prob){
  string2num = function(xstr) as.numeric(unlist(strsplit(xstr, " "))) 
  clust_names = names(clust_centers_w_prob$assign)
  n_clust = length(clust_names)
  oldpar = par(mar = (c(5,2,4,2)+.1))
  # hardcoded range, avoid?
  plot(x=slopes, y=rep(n_clust+1, length(slopes)), xlim =c(0,2),
       pch = 1, ylim = c(0,n_clust-1), yaxt = "n", ylab = NA,
       frame = F, xlab = "Fitted Slope", main = "Most Likely Clusters")
  axis(2, at = 0:n_clust, label = n_clust:0, las=1, lwd=0, lwd.ticks = 1)
  axis(4, at = 0:n_clust-1, label =paste("(", clust_centers_w_prob$assign[n_clust:0+1],")", sep=""),
       las=1, lwd = 0, line = -1)
  for(clust_id in 1: length(clust_names)){
    clust_assign = string2num(clust_names[clust_id])
    for( cid in unique(clust_assign)){
      points(x=slopes, y=rep(n_clust-clust_id, length(slopes)), xlim =c(min(slopes), 2),
             pch = ".", ylim = c(0,n_clust+1), yaxt = "n", ylab = NA)
      active_x = slopes[which(clust_assign==cid)]
      x_lims = range(active_x)+c(-.01, +.01)
      y_lims = rep(n_clust-clust_id, 2)+c(-.2, .2)
      rect(xleft =x_lims[1], xright = x_lims[2], 
           ybottom = y_lims[1], ytop = y_lims[2])
      #lines(x = c(range(active_x)), y=rep(n_clust-clust_id,2), lwd = 2)
    }
  }
  par(oldpar)
}


# Violin Plot for Scores ####
plot_scores = function(score_df,bootstrap_calc_list){
  plot_df = score_df
  names(plot_df) = c("Mix Desc", rep(names(bootstrap_calc_list), 3))
  n_calcs = length(bootstrap_calc_list)
  
  subdf = plot_df#[subindx,]
  mdf = rbind(cbind(melt(subdf[,grep("CRPS", names(score_df) )]), "Score" = "CRPS"),
              cbind(melt(subdf[,grep("MSE", names(score_df) )]), "Score" = "MSE"),
              cbind(melt(subdf[,grep("LLH", names(score_df))]), "Score" = "LLH"))
  names(mdf) = c("Method", "Value", "Score")
  
  ggplot(mdf, aes(x = Method, y=Value))+
    geom_boxplot(draw_quantiles = c(.25, .5, .75))+
    facet_wrap(vars(Score), scales = "free_y")+
    scale_y_log10()+
    labs(title = "Score Summary")
  
}

# Isobole plots ####
#create_isobole = function(a_vals = c(1,2), slope_vals = c(.5, 1, 1.5), use_3d = F){
if(F){  
  isobole_df = data.frame()
  for(a_val in c(-1, 1, 3, 5)){
    for(slope_val in c(.5, 1, 1.5)){
      # create a parameter input matrix 
      param_matrix = as.matrix(cbind("a" = c(a_val,5), 
                                     "b" = c(.8,1.5), 
                                     "c" = c(slope_val,1),
                                     "max_R" = 5,
                                     "d" = 0))
      # specify both chems in cluster 1 of 1
      cluster_assign_vec = c(1,1)
      # create a calculator to predict response given concentration
      mix_pred = mix_function_generator(param_matrix, cluster_assign_vec)
      # create matrix of concentration values
      c_val = seq(.1, 5, by=.1)
      n_pix = length(c_val)
      c_mat = expand.grid(c_val,c_val)
      # apply calculator to concentrations
      pred_effect = apply(c_mat, MARGIN=1, FUN = function(x) mix_pred(x) )
      mix_df = data.frame(cbind(c_mat, pred_effect))
      names(mix_df) = c("C1", "C2", "effect")
      isobole_df = rbind(isobole_df, cbind(mix_df, "a"=a_val, "slope" = slope_val))
    }
  }
  # 2d contour plot
  #pdf(file='Isoboles_RGCA_negative.pdf',width=8, height = 6)
  # ggplot(isobole_df, aes(x=C1, y=C2, z=effect)) + 
  #   geom_contour(bins = 25) + 
  #   facet_grid( rows = vars(slope),cols = vars(a))
  #dev.off()
  
  # 3d plot with 2d contour plot on floor
  #mix_df$effect[mix_df$effect< -1000]=NA
  plot_ly(x = (c_val),y=(c_val),z=matrix(mix_df$effect, nrow = n_pix)) %>% 
    add_surface(contours= list( z = list(
      show=TRUE,
      usecolormap=TRUE,
      highlightcolor="#ff0000",
      project=list(z=TRUE)
    ), zlim = c(0,1)
    )) %>% layout(scene = list(
      camera=list(eye = list(x=2, y=.5, z=.7))))
  #%>% add_surface(z = matrix(pred_regular, nrow = n_pix), opacity = .4) 
  #%>% add_surface(z = matrix(pred_agonist, nrow = n_pix), opacity = .4)
  #beepr::beep()
  
  #plot_ly(z=matrix(mix_df$effect, nrow = n_pix), type = "contour")
}


plot_individual_response = function(replicate_sets, re_par_summary,curve_fits, RE_curve_fits, Cx, y_i){
  replicate_counts = rep(1:18, unlist(lapply(replicate_sets, length)))
  # get MSE
  DRM_SSE = 0
  REM_SSE = 0
  # pdf(width=10, height = 8, onefile = F)
  old_par = par(mfrow=c(3, 3))
  for(ri in 1:length(replicate_sets)){
    u_RE_vals = re_par_summary$u_RE_params[which(replicate_counts==ri)]
    v_RE_vals = re_par_summary$v_RE_params[which(replicate_counts==ri)]
    
    rep_idx = replicate_sets[[ri]]
    plot(Cx[rep_idx,],y_i[rep_idx,], log="x", pch=1:length(rep_idx), 
         xlab = "Concentration", ylab = "Response", main = chem_map[ri])
    lines(Cx[rep_idx[1],], hill_function(curve_fits[ri,1],curve_fits[ri,2],curve_fits[ri,3],Cx[rep_idx[1],]),col=2)
    DRM_SSE = DRM_SSE + sum((y_i[rep_idx,] - hill_function(curve_fits[ri,1],curve_fits[ri,2],curve_fits[ri,3],Cx[rep_idx[1],]))^2, na.rm = T)
    for(re_idx in 1:length(u_RE_vals)){
      U_RE = u_RE_vals[re_idx]; V_RE = v_RE_vals[re_idx]
      RE_fit_pred = V_RE+hill_function(RE_curve_fits[ri,1]+U_RE,RE_curve_fits[ri,2],RE_curve_fits[ri,3],Cx[rep_idx[re_idx],])
      lines(Cx[rep_idx[re_idx],],RE_fit_pred ,col=3)
      REM_SSE = REM_SSE + sum((y_i[rep_idx[re_idx],] - RE_fit_pred)^2, na.rm = T)
    }
  }
  #dev.off()
  total_points = sum(!is.na(Cx))# length(unlist(replicate_sets))*ncol(Cx)
  DRM_SSE/total_points
  REM_SSE/total_points
  # return plot settings
  par(old_par)
  
}



# Hill plots ####

# creates the Hill function example with reflected segments used in the manuscript
plot_reflection_illustration = function(){
  require(scales);require(ggplot2)
  a = .9; b =.5; slope_c = 1.5
  #max_R = 2; d=0
  # max_R allows for variable definition of partial agonist
  # rather than assuming partial means < 1 response
  # a: sill parameter
  # b:  EC50 or Disassociation constant
  # c: slope
  # d: minimum effect
  hilly_inverse_test = function(y, c) hill_invs_factry(a,b,c)(y)
  #y_test_old = seq(abs(a)+.02, max_R-.01, by=.001)
  y_test = seq(a+.04, 2*a-.01, by=sign(a)*.01)
  y_test_bey = seq(2*a-.01, 6*a-.01, by=sign(a)*.01)
  y_test_neg = seq(0, -4*a-.01, by=-sign(a)*.01)
  x_test = seq(0, 5, by=.1)#1*10^(seq(-8, .5, by=.1)) 
  
  
  lgnd_brks = c(" GCA", "RGCA", "Extension", "Intermediate", "Reflection", "Ext+Reflect")
  
  pt1 = cbind(x_test, a/(1+(b/x_test)^(slope_c)), 1, 1, lgnd_brks[2])
  pt1_ghost = cbind(-x_test - 2*b, a/(1+(b/x_test)^(slope_c)), 5, 1, lgnd_brks[4])
  pt2 = cbind(sapply(y_test,function(x) hilly_inverse_test(x,slope_c)),
              y_test, 2, 1, lgnd_brks[5])
  pt3 = cbind(sapply(y_test_neg,function(x) hilly_inverse_test(x,slope_c)),
              y_test_neg, 3, 1, lgnd_brks[3])
  pt4 = cbind(sapply(y_test_bey,function(x) hilly_inverse_test(x,slope_c)),
              y_test_bey, 4, 1,lgnd_brks[6])
  x_test = seq(-5,5, by=.1)
  pt_gca = cbind(x_test, a/(1+(b/x_test)^(1)), 0, 2, lgnd_brks[1])
  pt_gca[which(pt_gca[,2]=="Inf"),2] = NA
  
  hill_df = data.frame(rbind(pt_gca,pt1, pt2, pt3, pt4,pt1_ghost))
  names(hill_df) = c("x", "y", "subcurve", "line", "Segment")
  hill_df$x=as.numeric(hill_df$x)
  hill_df$y=as.numeric(hill_df$y)
  
  rect_region = data.frame(
    xmin = c(0, -b),
    xmax = c(5, 0),
    ymin = c(0,-5),
    ymax= c(a, 0),
    z = c(1, 2)
  )
  
  okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  basic_col = c("black", "green","red", "gray", "cyan", "purple")
  ggplot(data = hill_df, aes(x=x,y=y, color = Segment))+ 
    geom_rect(data = rect_region, inherit.aes = FALSE, 
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = as.factor(z)),
              alpha = .2, show.legend = F)+
    geom_line(linewidth = 1, lineend = "round", linetype = as.numeric(hill_df$line))+
    geom_hline(yintercept = a, linetype = 2,color = "gray")+
    annotate("text", x=1, y=1, label= "Sill", size = 3, col="blue")+
    geom_vline(xintercept = -b, linetype = 2,color = "gray")+
    annotate("text", x=-.55, y=-1, label= "-EC50", angle = 90, size = 3, col="blue")+
    #   geom_abline(slope = -a/b)+
    #    geom_abline(slope = b/a, intercept = -(b/a*(-b)-a))+
    theme_bw()+
    theme(legend.position = "bottom",legend.box="vertical", legend.margin=margin())+#legend.justification=c(1,0), legend.position=c(1, .5))+#
    scale_colour_manual(values = basic_col, breaks = lgnd_brks)+
    xlab("Concentration (uM)")+
    coord_cartesian(xlim=c(-3, 2), ylim=c(-1.5, 3))+
    guides(colour = guide_legend(nrow = 1, override.aes =list(linetype = c(3, 1, 1,1, 1,1))))+
    ylab("Response")+labs(title = "Extending GCA with Reflection")#, subtitle = "")
  
  
  #pdf("RGCA_symmetry_detail.pdf", width = 8 ,height = 5)
  #png("RGCA_symmetry_full_neg.png",width = 8 ,height = 5, units = "in",res = 200)
  
}
#pdf("RGCA_symmetry_full.pdf", width =8,height = 4)
# plot_reflection_illustration()
# dev.off()

# not used, but plots the inverse of the Hill function 
plot_inverse_illustration = function(){
  require(scales);require(ggplot2)
  a = .9; b =.5; slope_c = 1.5
  #max_R = 2; d=0
  # max_R allows for variable definition of partial agonist
  # rather than assuming partial means < 1 response
  # a: sill parameter
  # b:  EC50 or Disassociation constant
  # c: slope
  # d: minimum effect
  hilly_inverse_test = function(y, c) hill_invs_factry(a,b,c)(y)
  #y_test_old = seq(abs(a)+.02, max_R-.01, by=.001)
  y_test = seq(a+.04, 2*a-.01, by=sign(a)*.01)
  y_test_bey = seq(2*a-.01, 6*a-.01, by=sign(a)*.01)
  y_test_neg = seq(0, -4*a-.01, by=-sign(a)*.01)
  x_test = seq(0, 5, by=.1)#1*10^(seq(-8, .5, by=.1)) 
  
  
  lgnd_brks = c(" GCA", "RGCA", "Extension", "Intermediate", "Reflection", "Ext+Reflect")
  
  pt1 = cbind(a/(1+(b/x_test)^(slope_c)),x_test, 1, 1, lgnd_brks[2])
  pt1_ghost = cbind(-x_test - 2*b, a/(1+(b/x_test)^(slope_c)), 5, 1, lgnd_brks[4])
  pt2 = cbind(y_test,
              sapply(y_test,function(x) hilly_inverse_test(x,slope_c)),
              2, 1, lgnd_brks[5])
  pt3 = cbind(y_test_neg,
              sapply(y_test_neg,function(x) hilly_inverse_test(x,slope_c)),
              3, 1, lgnd_brks[3])
  pt4 = cbind(y_test_bey,sapply(y_test_bey,function(x) hilly_inverse_test(x,slope_c)),
              4, 1,lgnd_brks[6])
  x_test = seq(-5,5, by=.1)
  #pt_gca = cbind(x_test, a/(1+(b/x_test)^(1)), 0, 2, lgnd_brks[1])
  pt_gca = cbind(a/(1+(b/x_test)^(1)),x_test,  0, 2, lgnd_brks[1])
  pt_gca[which(pt_gca[,1]=="Inf"),2] = NA
  
  hill_df = data.frame(rbind(pt_gca,pt1, pt2, pt3, pt4,pt1_ghost))
  names(hill_df) = c("x", "y", "subcurve", "line", "Segment")
  hill_df$x=as.numeric(hill_df$x)
  hill_df$y=as.numeric(hill_df$y)
  
  rect_region = data.frame(
    xmin = c(0, -b),
    xmax = c(5, 0),
    ymin = c(0,-5),
    ymax= c(a, 0),
    z = c(1, 2)
  )
  
  okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  basic_col = c("black", "green","red", "gray", "cyan", "purple")
  ggplot(data = hill_df, aes(x=x,y=y, color = Segment))+ 
    geom_line(linewidth = 1, lineend = "round", linetype = as.numeric(hill_df$line))+
    #geom_hline(yintercept = a, linetype = 2,color = "gray")+
    #annotate("text", x=1, y=1, label= "Sill", size = 3, col="blue")+
    #geom_vline(xintercept = -b, linetype = 2,color = "gray")+
    #annotate("text", x=-.55, y=-1, label= "-EC50", angle = 90, size = 3, col="blue")+
    #   geom_abline(slope = -a/b)+
    #    geom_abline(slope = b/a, intercept = -(b/a*(-b)-a))+
    theme_bw()+
    theme(legend.position = "bottom",legend.box="vertical", legend.margin=margin())+#legend.justification=c(1,0), legend.position=c(1, .5))+#
    scale_colour_manual(values = basic_col, breaks = lgnd_brks)+
    xlab("Concentration (uM)")+
    coord_cartesian(xlim=c(-3, 4), ylim=c(-3, 4))+
    guides(colour = guide_legend(nrow = 1, override.aes =list(linetype = c(3, 1, 1,1, 1,1))))+
    ylab("Response")+labs(title = "Extending GCA with Reflection, Inverse")#, subtitle = "")
  
  
  #pdf("RGCA_symmetry_detail.pdf", width = 8 ,height = 5)
  #png("RGCA_symmetry_full_neg.png",width = 8 ,height = 5, units = "in",res = 200)
  
}



# No Effect Inclusion Test #####
# Include/exclude no effects testing:
# for a mixture with chemicals that have 0 effect, does their
# inclusion affect the mixture response?
compare_exclude_include = function(){
  
  plot_compare = function(idx_set_list, conc_factor=1, plt_title = NA){
    reptit = sapply(1:length(idx_set_list), FUN = function(x) length(idx_set_list[[x]]))
    EC_df = cbind(mix_df[unlist(idx_set_list),], 
                  "group" = rep(names(idx_set_list), reptit),
                  "multp" = rep(conc_factor, reptit))
    # calibrate actual dosage using the guide
    plot_df = cbind(EC_df[,c(2, 4, 6)], EC_df[,c(6+15,51, 52) ])
    for(idxoff in 1:14){
      sub_df = cbind(EC_df[,c(2, 4, 6+idxoff)], EC_df[,c(6+idxoff+15,51, 52) ])
      names(sub_df) = names(plot_df)
      plot_df = rbind(plot_df,  sub_df)
    }
    plot_df$ReplicateSet = as.factor(plot_df$ReplicateSet)
    plot_df$CONC0 = plot_df$CONC0 * plot_df$multp
    p = ggplot(plot_df, aes(x = CONC0, y = RESP0, group = group, color = group, shape =ReplicateSet)) + 
      geom_point(size = 3, stroke=1.5)+ scale_x_continuous(trans='log10')+
      xlab("Concentration (uM)")+
      ylab("Response")+
      labs(title = plt_title)
    return(p)
  }
  
  
  
  EC_AR_guide = as.numeric(mix_guide[which(mix_guide$CAS ==  "NOCAS_49497"), 4:21])
  EC_all_guide = as.numeric(mix_guide[which(mix_guide$CAS ==  "NOCAS_49538"), 4:21])
  EC_ER_guide = as.numeric(mix_guide[which(mix_guide$CAS ==  "NOCAS_49517"), 4:21])
  crat_1 = EC_AR_guide/EC_all_guide
  mean(crat_1[crat_1>0])
  crat_2 = EC_ER_guide/EC_all_guide
  mean(crat_2[crat_2>0])
  
  idx_ec_ar = which(mix_df$CAS == "NOCAS_49497") # EquiConc AR
  idx_ec_all = which(mix_df$CAS == "NOCAS_49538") # EquiConc all
  idx_ec_er = which(mix_df$CAS == "NOCAS_49517") # EquiConc ER
  idx_ep_ar2 = which(mix_df$CAS == "NOCAS_49561") # EquiConc AR alt
  idx_set_list = list("EC_all" = idx_ec_all, "EC_AR" = idx_ec_ar, "EC_ER"= idx_ec_er, "EC_AR2" = idx_ep_ar2)
  p1=plot_compare(idx_set_list, conc_factor =c(1,1,3,1), plt_title = "Equiconcentration Comparison")
  
  EP_all_guide = as.numeric(mix_guide[which(mix_guide$CAS ==  "NOCAS_49537"), 4:21])
  EP_AR_guide= as.numeric(mix_guide[which(mix_guide$CAS ==  "NOCAS_49496"), 4:21])
  EP_ER_guide = as.numeric(mix_guide[which(mix_guide$CAS ==  "NOCAS_49516"), 4:21])
  crat_1 = EP_AR_guide/EP_all_guide
  mean(crat_1[crat_1>0])
  crat_2 = EP_ER_guide/EP_all_guide
  mean(crat_2[crat_2>0])
  
  idx_ep_all = which(mix_df$CAS == "NOCAS_49537") # EquiPot all
  idx_ep_ar = which(mix_df$CAS == "NOCAS_49496") # EquiPot AR
  idx_ep_er = which(mix_df$CAS == "NOCAS_49516") # EquiPot ER
  
  EP_set_list = list("EP_all" = idx_ep_all, "EP_AR"=idx_ep_ar, "EP_ER"=idx_ep_er)
  p2=plot_compare(EP_set_list, conc_factor =c(1,4,2), plt_title = "Equipotent Comparison")
  
  
  plot_grid(p1,p2)
  #EC_all and EC_AR may be mixed up!
}

# Sim Study ####
sim_study_boxplot = function(record_MSE, record_CRPS){
  # Plot results as boxplot
  record = data.frame(rbind(cbind(record_MSE, "Score" = "MSE"), cbind(record_CRPS,"Score" = "CRPS")))
  names(record) = c("RGCA+DP", "GCA", "IA", "Mix_Size", "Score")
  record_df = melt(record, id.vars = c("Score", "Mix_Size"))
  names(record_df) = c("Score_Method", "Mix_Size","Pred_Method", "Score")
  ggplot(record_df, aes(x = as.factor(Mix_Size), y = Score)) +
    geom_boxplot(aes(colour = Pred_Method))+
    facet_wrap(vars(Score_Method), scales= "free" )+
    theme_bw()+
    xlab("Number of Chemicals in Mixture")+
    guides(colour=guide_legend(title="Method"))
  
}