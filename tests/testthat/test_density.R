context("Density calculation")

test_that("density_calc output as expected", {
  expect_equal(density_calc_expt(), density_calc(density_calc_test(),
                                                radius = sqrt(1/pi)))
})

test_that("density_summary output as expected", {
  expect_equal(density_summ_expt(), density_summ_test(fake_map(), stand = "A"))
})

test_that("density_all_stands output as expected", {
  expect_equal(density_all_expt(), density_all_test(fake_map()))
})

test_that("density_specific output as expected", {
  expect_equal(density_spec_expt(), density_spec_test(fake_map(), stand = "A"))
})
