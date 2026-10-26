#fast regression tests: the same demo models as tests/full/test-fits.R, but fitted on coarse masks (see
#fit_demo_coarse() in helper-fits.R), so the expected values differ from the ones there

expected_fast = list(
  bearing_dist_hn = list(
    linked = c(`g0_link` = 36.04365339, `sigma_link` = 1.637706767, `kappa_link` = 2.543717688, `alpha_link` = 1.307177035, `D_link` = 7.786664818),
    fitted = c(`g0` = 1, `sigma` = 5.143361052, `kappa` = 12.72689776, `alpha` = 3.695726067, `D` = 2408.272148),
    se = c(`g0_link` = 0, `sigma_link` = 0.03320028128, `kappa_link` = 0.12299169, `alpha_link` = 0.09638686666, `D_link` = 0.1042061795)),
  ihd_ext = list(
    linked = c(`g0_link` = 36.04365339, `sigma.(Intercept)_link` = 2.014238054, `sigma.brandsony_link` = -0.2688539972, `D.(Intercept)_link` = 5.791007404, `D.noise_link` = 0.1208865741),
    fitted = c(`g0` = 1, `sigma.(Intercept)` = 2.014238054, `sigma.brandsony` = -0.2688539972, `D.(Intercept)` = 5.791007404, `D.noise` = 0.1208865741),
    se = c(`g0_link` = 0, `sigma.(Intercept)_link` = 0.1346702388, `sigma.brandsony_link` = 0.06917416247, `D.(Intercept)_link` = 2.061305657, `D.noise_link` = 0.2103649803)),
  mul_ses_ext = list(
    linked = c(`g0.(Intercept)_link` = 1.243709083, `g0.weathersunny_link` = -0.8880297397, `sigma.(Intercept)_link` = 1.104038442, `sigma.brandsony_link` = 0.1678180405, `kappa_link` = 2.536984614, `alpha_link` = 0.9034012618, `D_link` = 7.884606598),
    fitted = c(`g0.(Intercept)` = 1.243709083, `g0.weathersunny` = -0.8880297397, `sigma.(Intercept)` = 1.104038442, `sigma.brandsony` = 0.1678180405, `kappa` = 12.64149445, `alpha` = 2.467983108, `D` = 2656.079914),
    se = c(`g0.(Intercept)_link` = 1.249053889, `g0.weathersunny_link` = 1.05225023, `sigma.(Intercept)_link` = 0.08145683814, `sigma.brandsony_link` = 0.07908663065, `kappa_link` = 0.3939801968, `alpha_link` = 0.18172571, `D_link` = 0.182429707)),
  hhn_cue = list(
    linked = c(`sigma_link` = 1.407863893, `lambda0_link` = 1.518207805, `D_link` = 6.199604056),
    fitted = c(`sigma` = 4.087215343, `lambda0` = 4.564038215, `D` = 492.5539787),
    se = c(`sigma_link` = 0.07805275914, `lambda0_link` = 0.1972011764, `D_link` = 0.1369068919)),
  ss_toa = list(
    linked = c(`b0.ss_link` = 4.503695747, `b1.ss_link` = 1.329525025, `sigma.ss_link` = 2.24474785, `sigma.toa_link` = -5.861862031, `D_link` = 7.744322514),
    fitted = c(`b0.ss` = 90.35042734, `b1.ss` = 3.779247911, `sigma.ss` = 9.438035453, `sigma.toa` = 0.002845939518, `D` = 2308.429065),
    se = c(`b0.ss_link` = 0.01474176824, `b1.ss_link` = 0.05342008424, `sigma.ss_link` = 0.05137365757, `sigma.toa_link` = 0.07677098329, `D_link` = 0.1046955126)),
  ind_toa_hhn = list(
    linked = c(`sigma_link` = 0.7958189158, `lambda0_link` = 1.800216065, `sigma.toa_link` = -6.549855281, `D.(Intercept)_link` = 10.60940501, `D.noise_link` = -0.4909001126, `mu_link` = 2.137970691),
    fitted = c(`sigma` = 2.21625518, `lambda0` = 6.050954724, `sigma.toa` = 0.001430322579, `D.(Intercept)` = 10.60940501, `D.noise` = -0.4909001126, `mu` = 8.482207136),
    se = c(`sigma_link` = 0.0812713723, `lambda0_link` = 0.2672819056, `sigma.toa_link` = 0.09869009826, `D.(Intercept)_link` = 5.035568297, `D.noise_link` = 0.4715096158, `mu_link` = 0.1530350522)),
  ind_ss_log = list(
    linked = c(`b0.ss_link` = 1.67482626, `b1.ss_link` = -2.470503553, `sigma.ss_link` = 2.356111826, `D_link` = 4.564630112, `mu_link` = 2.165748184),
    fitted = c(`b0.ss` = 5.337867665, `b1.ss` = 0.08454227676, `sigma.ss` = 10.54985195, `D` = 96.02706821, `mu` = 8.721124481),
    se = c(`b0.ss_link` = 0.001713788785, `b1.ss_link` = 0.01761343554, `sigma.ss_link` = 0.02468066938, `D_link` = 0.1743767421, `mu_link` = 0.06392061258)),
  ind_ss_sp = list(
    linked = c(`b0.ss_link` = 4.628315967, `b1.ss_link` = -1.664330401, `sigma.ss_link` = 2.30526273, `D_link` = 4.240193206, `mu_link` = 2.247590948),
    fitted = c(`b0.ss` = 102.3415723, `b1.ss` = 0.1893173822, `sigma.ss` = 10.02681225, `D` = 69.42126312, `mu` = 9.464906897),
    se = c(`b0.ss_link` = 0.02068886311, `b1.ss_link` = 0.6090011913, `sigma.ss_link` = 0.0279840291, `D_link` = 0.2000002474, `mu_link` = 0.06543279516))
)

expect_fit_matches = function(fit, expected, tolerance = 1e-4){
  expect_equal(coef(fit, types = 'linked'), expected$linked, tolerance = tolerance)
  expect_equal(coef(fit, types = 'fitted'), expected$fitted, tolerance = tolerance)
  expect_equal(stdEr(fit, types = 'linked'), expected$se, tolerance = tolerance)
}

for(data_name in names(expected_fast)){
  test_that(paste0("fast fit on coarse mask -- ", data_name), {
    fit = fit_demo_coarse(data_name)
    expect_fit_matches(fit, expected_fast[[data_name]])
  })
}

test_that("gradient-free fit on coarse mask agrees with the gradient-based one", {
  fit = fit_demo_coarse('ihd_ext', gradient_free = TRUE)
  expect_fit_matches(fit, expected_fast[['ihd_ext']])
})

test_that("confidence intervals are consistent with estimates and std errors", {
  fit = fit_demo_coarse('hhn_cue')
  ci = confint(fit, types = 'linked', level = 0.9)
  est = coef(fit, types = 'linked')[rownames(ci)]
  se = stdEr(fit, types = 'linked')[rownames(ci)]
  expect_equal(unname(ci[, 1]), unname(est + stats::qnorm(0.05) * se), tolerance = 1e-6)
  expect_equal(unname(ci[, 2]), unname(est + stats::qnorm(0.95) * se), tolerance = 1e-6)
  expect_equal(colnames(ci), c('5 %', '95 %'))
})
