library(dplyr)
library(tidyr)

# functions to deal with Tox21 data:
# function to convert named array to vector
mk_vec <- function(arow) unlist(array(arow))
# functions to take dot or matrix products when nulls present:
# vector prod
vprod_NA <- function(a, b) {
  mat_out <- (t(na.omit(as.matrix(a))) %*% na.omit(as.matrix(b)))
  if (dim(a)[2] > 1 || dim(b)[2] > 1) {
    return(mat_out)
  }
  return(mat_out[1])
}
# quadratic product
qprod_NA <- function(a, b, c) {
  return((t(na.omit(as.matrix(a))) %*% diag(na.omit(b)) %*% na.omit(as.matrix(c)))[1])
}
# outer product
oprod_NA <- function(a, b) {
  return(outer(na.omit(a), na.omit(b)))
}



# read in data
# notes about input data
# AR-luc is best in terms of not having synergy/antagonism and no dose issues
# AR-bla seems to have a lot of antagonism
# ER-luc has a dose issue with zearalenone, throwing off a lot of estimates
# ER-bla (same issue as ER-luc?)
#


#' Read preprocessed Tox21 Input data
#'
#' @param input_file processed data sets from Tox21 with mixtures, included
#' with package
#'
#' @return
#' @export
#'
#' @examples
read_prepared_Tox21_data <- function(input_file = "AR-luc.txt") {
  df <- read.delim(paste("Input/", input_file, sep = ""))
  # head(df)
  # pull mix data for computing concentrations of mixture components
  mix_df <<- df %>% filter(substr(CAS, 1, 5) == "NOCAS")
  # separate pure data for fitting tox model
  pure_df <<- df %>% filter(substr(CAS, 1, 5) != "NOCAS") # %>% filter(Curve.Class!=4)
  # index of replicates
  cas_rep_idx <- c(2, 4)
  conc_idx_T21_matrix <- 6:20 # concentration values (x)
  conc_colnames <- colnames(pure_df)[conc_idx_T21_matrix]
  resp_idx <- 21:35 # response values (y)
  resp_colnames <- colnames(pure_df)[resp_idx]
  mask_idx <- 36:50 # masked values due to cytotoxicity, etc
  mask_colnames <- colnames(pure_df)[mask_idx]
  # correct for duplicated replicates:  add 3
  dup_replicates <- duplicated(pure_df[, cas_rep_idx])
  pure_df$ReplicateSet[dup_replicates] <- pure_df$ReplicateSet[dup_replicates] + 3

  Cx_T21 <- pure_df[, c(cas_rep_idx, conc_idx_T21_matrix)]
  y_i_T21 <- pure_df[, c(cas_rep_idx, resp_idx)]
  mask_T21 <- cbind(pure_df[, cas_rep_idx], pure_df[, mask_idx] == T)

  # simple replicate handling:  rowbind
  # y_i_T21reps = y_i_T21 %>%
  #   pivot_wider(names_from = "ReplicateSet", values_from = all_of(resp_colnames) )
  #
  # Cx_T21reps = Cx_T21 %>%
  #   pivot_wider(names_from = "ReplicateSet", values_from = all_of(conc_colnames) )
  #
  # mask_reps = mask_T21 %>%
  #   pivot_wider(names_from = "ReplicateSet", values_from = all_of(mask_colnames) )


  # better replicate handling:  merge CAS and replicate ID and specify random effect
  n_chems <<- length(unique(pure_df$CAS))
  rep_table <- table(pure_df$CAS)
  n_replicates_per_chem <<- unlist(array(rep_table))
  chem_id <<- names(rep_table)

  replicate_sets <<- c()
  pure_unique_CAS <- chem_id # unique(y_i_T21$CAS)


  for (chem in pure_unique_CAS) {
    rep_idx <- which(y_i_T21$CAS == chem)
    replicate_sets <<- c(replicate_sets, list(rep_idx))
  }

  mask_layer <- as.matrix(mask_T21[, 3:ncol(mask_T21)])
  y_i <- y_i_T21[, 3:ncol(y_i_T21)]
  y_i[mask_layer] <- NA
  y_i <<- as.matrix(unname(y_i))
  Cx <- Cx_T21[, 3:ncol(y_i_T21)]
  Cx[mask_layer] <- NA
  Cx <<- as.matrix(unname(Cx))
  # scale control from 100 to 1
  # y_i = y_i / 100




  # Pull out mixtures:  need mixdf for effect,
  mix_guide <<- readxl::read_xls("Input/AllMixtureComponentsARER.xls")
  # need chemical order, apply to concentration

  CAS_nums <<- colnames(mix_guide)[4:21]
  c_cm <- mix_guide[, 4:21]
  mix_conc_df <<- mix_df[, conc_idx_T21_matrix]
  # extrapolate beyond?
  extrap_df <- mix_conc_df * 1e5
  names(extrap_df) <- paste(names(extrap_df), "xt", sep = "")
  mix_conc_df <<- cbind(mix_conc_df, extrap_df)
  rel_mix_conc_df <- sweep(mix_conc_df,
    MARGIN = 1,
    FUN = "/",
    STATS = mix_conc_df$CONC14
  )

  # create a function to access the mixture data
  get_conc_matrix <<- function(row_idx) {
    # for specific mix, use  cm*crel*.01
    guide_idx <- which(mix_guide$CAS == mix_df[row_idx, ]$CAS)
    n_tot_chems <- length(c_cm[guide_idx, ])
    # n_active_chems = n_chems = chems with function cat!=4
    # build matrix:  rows = conc, col = chemical
    conc_mix_matrix <- matrix(0,
      nrow = ncol(rel_mix_conc_df),
      ncol = n_tot_chems
    )
    # iterate across chemicals to build conc matrix
    for (chem_conc_idx in 1:n_tot_chems) {
      chem_conc_ra <- as.numeric(c_cm[guide_idx, chem_conc_idx])
      # multiply: conc of guide * relative dose  * 0.01
      chem_mix_conc <- chem_conc_ra * rel_mix_conc_df[row_idx, ] * 1e-6 * 1e-2
      # save as column
      conc_mix_matrix[, chem_conc_idx] <- unlist(chem_mix_conc)
    }

    return(conc_mix_matrix)
  }


  # reference to identify curves
  AR_agonist_rows <<- c(6, 7, 9, 10, 11, 12, 13, 18)
  ER_agonist_rows <<- dplyr::setdiff(1:18, AR_agonist_rows)
  # pure unique ordering for plotting ####
  relevant_guide_CASs <<- which(CAS_nums %in% pure_unique_CAS)
  # reorder the concentration from the guide to match the read-in data
  pure_ordering <<- sapply(pure_unique_CAS, FUN = function(x) which(x == CAS_nums[relevant_guide_CASs]))
  chem_map_plain <- sapply(names(pure_ordering), FUN = function(cas_num) pure_df$Sample.Name[which(pure_df$CAS == cas_num)[1]])
  chem_map <<- chem_map_plain
  chem_map[AR_agonist_rows] <<- paste("*", chem_map[AR_agonist_rows])
  names(pure_ordering) <<- chem_map
}


#' Generate a string based on the vector of mixture indices that 
#' satisfy the boolean input conditions
#'
#' @param drop_no_effect fix agonists for other assays to 0
#' @param use_graded_filter mixtures with the term "graded" in description
#' @param filter_graded_with_all_nonzero graded mixtures involving all chemicals
#' @param use_binary_filter mixtures with the term "binary" in description
#' @param ignore_no_effect_mix if all AR chemicals have 0 concentration
#' @param skip_zearalenone Data issues for ER data, skipped
#'
#' @return
#' @export
#'
#' @examples
generate_mixture_subset <- function(drop_no_effect = F,
                                    use_graded_filter = T,
                                    filter_graded_with_all_nonzero = F,
                                    use_binary_filter = F,
                                    ignore_no_effect_mix = F,
                                    skip_zearalenone = F) {
  # collect which mixtures have desired descriptions
  valid_idx <- c()
  for (mix_idx in 1:nrow(mix_df)) {
    # skip replicates
    skip_plot <- ifelse(mix_df$ReplicateSet[mix_idx] > 1, T, F)
    if (skip_plot) next

    # check for graded mixtures and equiconcentration
    CAS_desc <- mix_guide$Description[which(mix_guide$CAS == mix_df$CAS[mix_idx])]
    is_graded <- any(grep("Graded", CAS_desc, ignore.case = TRUE))
    is_equicon <- any(grep("EQUICONCENTRATION", CAS_desc, ignore.case = TRUE))
    is_binary <- any(grep("binary", CAS_desc, ignore.case = TRUE))
    contains_zearalenone <- mix_guide$`17924-92-4`[which(mix_guide$CAS == mix_df$CAS[mix_idx])] > 0
    chem_conc_matr <- get_conc_matrix(mix_idx)

    if (skip_zearalenone) {
      if (contains_zearalenone) next
    }

    if (use_graded_filter) {
      if (!(is_graded || is_equicon)) next
      # skip equiconcentration if does not include all chems
      if (is_equicon && any(chem_conc_matr[1, ] == 0)) next
      # skip graded if it does not include all chems
      if (filter_graded_with_all_nonzero && is_graded && any(chem_conc_matr[1, ] == 0)) next
    }
    if (ignore_no_effect_mix) {
      if (all(chem_conc_matr[1, AR_agonist_rows] == 0)) next
    }
    if (drop_no_effect) chem_conc_matr[, ER_agonist_rows] <- 0
    if (use_binary_filter) {
      if (!is_binary) next
    }

    valid_idx <- c(valid_idx, mix_idx)
  }
  valid_idx
  return(paste(valid_idx, collapse = ","))
}
