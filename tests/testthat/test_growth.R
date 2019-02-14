context("Annual growth calculations")

test_that("growth_summary output as expected", {
  expect_equal(growth_summ_exp(), as.data.frame(growth_summary(fake_growth())))
})

test_that("detailed_growth output as expected", {
  expect_equal(det_growth_exp(), as.data.frame(detailed_growth(fake_growth())))
})

test_that("defined_period_annual_growth output as expected", {
  expect_equal(def_growth_exp(), 
               as.data.frame(defined_period_annual_growth(fake_growth(),
                                                          2005, 2010)))
})