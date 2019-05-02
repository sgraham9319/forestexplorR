context("Density calculation")

test_that("density_calc output as expected", {
  expect_equal(density_calc_exp(), density_calc(density_calc_test(),
                                                radius = sqrt(1/pi)))
})

test_that("density_summary output as expected", {
  expect_equal(density_summ_exp(), density_summ_test(fake_map(), stand = "A"))
})

test_that("density_all_stands output as expected", {
  expect_equal(density_all_exp(), density_all_test(fake_map()))
})

test_that("density_specific output as expected", {
  expect_equal(density_spec_exp(), density_spec_test(fake_map(), stand = "A"))
})
