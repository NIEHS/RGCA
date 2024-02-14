test_that("rgca hill inverse works", {
  a <- 2
  b <- 3
  c <- 0.5
  y <- -1
  # case 1, y extended to small negative conc, invert slope
  expect_equal(hill_invs_factry(a, b, c)(y),
               -b / (1 + (-a / y)^(1 / c)))
  # negate a and y
  a <- -a
  y <- -y
  expect_equal(hill_invs_factry(a, b, c)(y),
               -b / (1 + (-a / y)^(1 / c)))
  # case 2, standard inverse, y<a
  a <- 2
  y <- 1
  expect_equal(hill_invs_factry(a, b, c)(y),
               b / (a / y - 1)^(1 / c))
  # negate a and y
  a <- -a
  y <- -y
  expect_equal(hill_invs_factry(a, b, c)(y),
               b / (a / y - 1)^(1 / c))
  # case 3, reflected part of the standard inverse, y < 2a
  a <- 2
  y <- 3.5
  expect_equal(hill_invs_factry(a, b, c)(y),
               -2 * b - b / (a / (2 * a - y) - 1)^(1 / c))
  # negate a and y
  a <- -a
  y <- -y
  expect_equal(hill_invs_factry(a, b, c)(y),
               -2 * b - b / (a / (2 * a - y) - 1)^(1 / c))
  # case 4, reflection of extension, slope inverted
  a <- 2
  y <- 6
  expect_equal(hill_invs_factry(a, b, c)(y),
               -2 * b + b / (1 + (a / (-2 * a + y))^(1 / c)))
  y <- -y
  a <- -a
  expect_equal(hill_invs_factry(a, b, c)(y),
               -2 * b + b / (1 + (a / (-2 * a + y))^(1 / c)))
})
