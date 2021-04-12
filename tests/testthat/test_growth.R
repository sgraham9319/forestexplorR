context("Annual growth calculations")

test_that("growth_summary output as expected", {
  expect_equal(growth_summ_expt(), growth_summ_test(fake_growth()))
})

test_that("detailed_growth output as expected", {
  expect_equal(det_growth_expt(), det_growth_test(fake_growth()))
})

test_that("defined_period_annual_growth output as expected", {
  expect_equal(def_growth_expt(), def_growth_test(fake_growth()))
})