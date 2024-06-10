test_that("rgca hill inverse works", {
  a <- 2
  b <- 3
  c <- 0.5
  y <- -1
  # case 0
  expect_equal(hill_invs_factry(a, b, c)(0), 0)
  # case 1, y extended to small negative conc, invert slope
  expect_equal(
    hill_invs_factry(a, b, c)(y),
    -b / (1 + (-a / y)^(1 / c))
  )
  # negate a and y
  a <- -a
  y <- -y
  expect_equal(
    hill_invs_factry(a, b, c)(y),
    -b / (1 + (-a / y)^(1 / c))
  )
  # case 2, standard inverse, y<a
  a <- 2
  y <- 1
  expect_equal(
    hill_invs_factry(a, b, c)(y),
    b / (a / y - 1)^(1 / c)
  )
  # negate a and y
  a <- -a
  y <- -y
  expect_equal(
    hill_invs_factry(a, b, c)(y),
    b / (a / y - 1)^(1 / c)
  )
  # case 3, reflected part of the standard inverse, y < 2a
  a <- 2
  y <- 3.5
  expect_equal(
    hill_invs_factry(a, b, c)(y),
    -2 * b - b / (a / (2 * a - y) - 1)^(1 / c)
  )
  # negate a and y
  a <- -a
  y <- -y
  expect_equal(
    hill_invs_factry(a, b, c)(y),
    -2 * b - b / (a / (2 * a - y) - 1)^(1 / c)
  )
  # case 4, reflection of extension, slope inverted
  a <- 2
  y <- 6
  expect_equal(
    hill_invs_factry(a, b, c)(y),
    -2 * b + b / (1 + (a / (-2 * a + y))^(1 / c))
  )
  y <- -y
  a <- -a
  expect_equal(
    hill_invs_factry(a, b, c)(y),
    -2 * b + b / (1 + (a / (-2 * a + y))^(1 / c))
  )
})

test_that("mix_response_prediction_works", {
  # specify a set of chemicals
  sills <- c(6, 3, 4, 8)
  ec50_vec <- c(1.2, 2.1, 0.1, 4.1)
  slopes <- c(0.5, 1.1, 2.0, 1.2)
  # Rmax is used to scale IA across clusters, can copy sills
  param_matrix <- as.matrix(cbind(
    "a" = sills,
    "b" = ec50_vec,
    "c" = slopes,
    "max_R" = sills
  ))

  # create the inverse function list used in denominator of GCA
  hill_inverse_list <- apply(param_matrix,
    MARGIN = 1,
    function(x) {
      do.call(hill_invs_factry, as.list(x))
    }
  )
  # create the GCA function to optimize over
  GCA_function <- eff_response_opt(hill_inverse_list,
    conc_vec = c(1, 2, 3),
    synergy_const = 0,
    interval_sign = 1
  )
  # Check GCA equation solving, match saved state (snapshot June 2024)
  expect_snapshot(optimize(GCA_function,
    interval = c(-20, 10),
    tol = .Machine$double.eps
  ))
  # test the negative case as well
  GCA_function_neg <- eff_response_opt(hill_inverse_list,
    conc_vec = c(1, 2, 3),
    synergy_const = 0,
    interval_sign = -1
  )
  expect_snapshot(optimize(GCA_function_neg,
    interval = c(-20, 10),
    tol = .Machine$double.eps
  ))

  # specify a clustering
  clust_assign <- c(1, 1, 2, 1)
  # generate the mix response function
  mix_function <- mix_function_generator(param_matrix,
    clust_assign,
    get_counts = FALSE,
    scale_CA = FALSE,
    synergy_const = 0
  )
  # specify some mixture doses to test
  dose_range <- seq(0, 10, length.out = 3)
  # create a matrix of mixture doses
  dose_matrix <- rev(expand.grid(
    dose_range, dose_range,
    dose_range, dose_range
  ))
  # test that mixture doses give expected response (snapshot April 2024)
  expect_snapshot(apply(dose_matrix, MARGIN = 1, mix_function))
})


test_that("summary_stats_are_correct", {
  set.seed(100)
  # create fake mcmc chains
  fake_MCMC <- list()
  iters <- 300
  n_chems <- 10
  n_reps <- 3
  fake_MCMC$slope_record <- matrix(abs(rnorm(n_chems * iters)), nrow = iters)
  fake_MCMC$sill_mideffect_record <- matrix(abs(rnorm(2 * n_chems * iters)),
    nrow = iters
  )
  fake_MCMC$sigma <- matrix(abs(rnorm(n_chems * iters)), nrow = iters)
  fake_MCMC$tau <- matrix(abs(rnorm(iters)), nrow = iters)
  fake_MCMC$u_RE <- matrix(rnorm(n_reps * n_chems * iters), nrow = iters)
  fake_MCMC$v_RE <- matrix(rnorm(n_reps * n_chems * iters), nrow = iters)
  fake_MCMC$u_RE_sd <- matrix(abs(rnorm(n_chems * iters)), nrow = iters)
  fake_MCMC$v_RE_sd <- matrix(abs(rnorm(n_chems * iters)), nrow = iters)

  # the summary stats should correspond to the median of the thinned samples
  param_summary <- pull_summary_parameters(fake_MCMC)
  expect_snapshot(param_summary)

  # prepare a cluster assignment
  clust_name <- paste(rep(c(1, 2), 5), collapse = " ")
  cluster_weight <- 1
  names(cluster_weight) <- clust_name
  clust_list <- list("cluster_assign" = cluster_weight)
  
  # for the second parameter summary function, need replicate sets
  repl_fun <- function(idx) idx + ((1:n_reps) - 1) * 10
  replicate_sets <<- lapply(1:n_chems, repl_fun)
  # the pulled parameters should correspond to the thinned samples only
  param_list <- pull_parameters(fake_MCMC, input_replicates = replicate_sets)
  expect_snapshot(param_list)

  # test the mix calculator generator using the param list
  param_list <- c(param_list, clust_list)
  mix_calc_test <- create_mix_calc(idx = 1, par_list = param_list)
  # test that the prediction is a number
  expect_equal(typeof(mix_calc_test(c(1:n_chems))), "double")
})

test_that("random_cluster_reproducible_with_seed", {
  set.seed(100)
  expect_equal(random_assignment(5, 3), c(1, 2, 1, 2, 3))
})

test_that("predicting_many_mix_responses_produces_correct_dim", {
  n_dose <- 10
  n_bootstraps <- 20
  n_chems <- 3
  chem_conc_matrix <- matrix(rnorm(n_chems * n_dose),
    nrow = n_dose,
    ncol = n_chems
  )
  bootstrap_calc_list <- lapply(1:n_bootstraps, FUN = function(x) {
    return(function(...) {
      return(1)
    })
  })
  bootstrap_calc_list <- list(bootstrap_calc_list, bootstrap_calc_list)
  n_methods <- length(bootstrap_calc_list)
  predict_resp_matrix <- predict_mix_response_many(
    n_dose,
    chem_conc_matrix,
    bootstrap_calc_list
  )
  expect_type(predict_resp_matrix, "list")
  expect_equal(length(predict_resp_matrix), n_methods)
  expect_equal(dim(predict_resp_matrix[[1]]), c(n_bootstraps, n_dose))
})

test_that("MLE_curve_fit_returns_estimates", {
  doses <- 10
  replicate_sets <- list(c(1, 2), c(3, 4), c(5, 6), c(7, 8))
  repls <- 2
  chems <- length(replicate_sets)
  set.seed(1)
  y_i <- matrix(abs(rnorm(doses * chems * repls)),
    nrow = repls * chems,
    ncol = doses
  )
  Cx <- matrix(rep(1:10, each = repls * chems),
    nrow = repls * chems,
    ncol = doses
  )
  drc_fit <- get_mle_curve_fits(y_i, Cx, replicate_sets)
  # check that each chemical got estimates for 3 parameters
  expect_equal(dim(drc_fit), c(chems, 3))
})
