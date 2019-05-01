context("Annual growth calculations")

test_that("growth_summary output as expected", {
  expect_equal(growth_summ_exp(), growth_summ_test(fake_growth()))
})

test_that("detailed_growth output as expected", {
  expect_equal(det_growth_exp(), det_growth_test(fake_growth()))
})

test_that("defined_period_annual_growth output as expected", {
  expect_equal(def_growth_exp(), def_growth_test(fake_growth()))
})