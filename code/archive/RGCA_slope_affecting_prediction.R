


# TODO minimum working example where slope affects prediction; used browser before


cluster_assign_vec = c(1,1,1,2,2,2,3,3,3)

# generate calculator for 
sills =  runif(18, max = 40, min=10)
param_matrix = as.matrix(cbind("a" =sills, 
                                         "b" =re_par_summary$ec50_params, 
                                         "c" = rep(1, 18),
                                         "max_R" =max(sills),
                                         "d" = 0))
mixture_calculator = mix_function_generator(bootstrap_param_matrix, 
                                            cluster_assign_vec)


create_mix_calc(x, tot_par_list))

mix_idx = 10
chem_conc_matr = get_conc_matrix(mix_idx)
if(drop_no_effect) chem_conc_matr[,bad_rows] = 0
n_dose = nrow(chem_conc_matr)-10 # first 15 are real, second 15 are extrapolation.  Here, just extrapolating lesss
chem_conc_matr_reordered = chem_conc_matr[,pure_ordering]
# each row of conc_matrix is a mixture dose
curve_data = matrix(0, ncol = n_dose, nrow = n_bootstraps)
curve_data_GCA = matrix(0, ncol = n_dose, nrow = n_bootstraps)
curve_data_IA = matrix(0, ncol = n_dose, nrow = n_bootstraps)

  conc_val = chem_conc_matr_reordered[1,]#+1e-40
  #Rprof()
  print(sampled_mix_funs[[1]](conc_val))
  #Rprof(NULL)
  #summaryRprof("Rprof.out")
  print(sampled_mix_funs_IA[[1]](conc_val))
