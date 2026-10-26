## Fast regression tests: demo models fitted on coarse masks (see fit.demo.coarse() in helper-fits.R).
## The expected values are specific to these masks.

expected.fast <- list(
  bearing_dist_hn = list(
    linked = c(`g0_link` = 36.04365339, `sigma_link` = 1.664843974, `kappa_link` = 2.176781347, `alpha_link` = 1.194431231, `D_link` = 7.742785896),
    fitted = c(`g0` = 1, `sigma` = 5.284848609, `kappa` = 8.817878838, `alpha` = 3.301679345, `D` = 2304.884617),
    se.linked = c(`g0_link` = 0, `sigma_link` = 0.03332972297, `kappa_link` = 0.1190753377, `alpha_link` = 0.09763469669, `D_link` = 0.1045368032),
    se.fitted = c(`g0` = 0, `sigma` = 0.1761425401, `kappa` = 1.049991901, `alpha` = 0.3223584614, `D` = 240.9452695)),
  ihd_ext = list(
    linked = c(`g0_link` = 36.04365339, `sigma.(Intercept)_link` = 2.015201369, `sigma.brandsony_link` = -0.2685835651, `D.(Intercept)_link` = 5.802687553, `D.noise_link` = 0.1195906146),
    fitted = c(`g0` = 1, `sigma.(Intercept)` = 2.015201369, `sigma.brandsony` = -0.2685835651, `D.(Intercept)` = 5.802687553, `D.noise` = 0.1195906146),
    se.linked = c(`g0_link` = 0, `sigma.(Intercept)_link` = 0.1340566375, `sigma.brandsony_link` = 0.0689811146, `D.(Intercept)_link` = 2.052882277, `D.noise_link` = 0.2094313595),
    se.fitted = c(`g0` = 0, `sigma.(Intercept)` = 0.1340566375, `sigma.brandsony` = 0.0689811146, `D.(Intercept)` = 2.052882277, `D.noise` = 0.2094313595)),
  mul_ses_ext = list(
    linked = c(`g0.(Intercept)_link` = 1.25359777, `g0.weathersunny_link` = -0.8911804938, `sigma.(Intercept)_link` = 1.100837189, `sigma.brandsony_link` = 0.1752579176, `kappa_link` = 2.576412486, `alpha_link` = 0.9306624646, `D_link` = 7.881289518),
    fitted = c(`g0.(Intercept)` = 1.25359777, `g0.weathersunny` = -0.8911804938, `sigma.(Intercept)` = 1.100837189, `sigma.brandsony` = 0.1752579176, `kappa` = 13.14987806, `alpha` = 2.536188757, `D` = 2647.284081),
    se.linked = c(`g0.(Intercept)_link` = 1.249589952, `g0.weathersunny_link` = 1.056282017, `sigma.(Intercept)_link` = 0.08094759863, `sigma.brandsony_link` = 0.07862947605, `kappa_link` = 0.2832507518, `alpha_link` = 0.1819861236, `D_link` = 0.1821153843),
    se.fitted = c(`g0.(Intercept)` = 1.249589952, `g0.weathersunny` = 1.056282017, `sigma.(Intercept)` = 0.08094759863, `sigma.brandsony` = 0.07862947605, `kappa` = 3.724712848, `alpha` = 0.4615511606, `D` = 482.1111577)),
  hhn_cue = list(
    linked = c(`sigma_link` = 1.395335085, `lambda0_link` = 1.514337902, `D_link` = 6.218601322),
    fitted = c(`sigma` = 4.036326857, `lambda0` = 4.546409963, `D` = 502.0006038),
    se.linked = c(`sigma_link` = 0.08283691265, `lambda0_link` = 0.1918478291, `D_link` = 0.1412736612),
    se.fitted = c(`sigma` = 0.3343568553, `lambda0` = 0.8722188816, `D` = 70.91946323)),
  ss_toa = list(
    linked = c(`b0.ss_link` = 4.488967996, `b1.ss_link` = 1.341111752, `sigma.ss_link` = 2.250698464, `sigma.toa_link` = -5.777377271, `D_link` = 7.813259924),
    fitted = c(`b0.ss` = 89.02951962, `b1.ss` = 3.823291696, `sigma.ss` = 9.49436499, `sigma.toa` = 0.003096826905, `D` = 2473.179687),
    se.linked = c(`b0.ss_link` = 0.01422445165, `b1.ss_link` = 0.05394471267, `sigma.ss_link` = 0.05099054628, `sigma.toa_link` = 0.06886454058, `D_link` = 0.1034190325),
    se.fitted = c(`b0.ss` = 1.266396097, `b1.ss` = 0.206246372, `sigma.ss` = 0.4841228575, `sigma.toa` = 0.0002132615621, `D` = 255.7738505)),
  ind_toa_hhn = list(
    linked = c(`sigma_link` = 0.9007327179, `lambda0_link` = 1.584816037, `sigma.toa_link` = -6.211909958, `D.(Intercept)_link` = 11.1479626, `D.noise_link` = -0.5586893986, `mu_link` = 2.296091984),
    fitted = c(`sigma` = 2.461405967, `lambda0` = 4.878393854, `sigma.toa` = 0.002005403567, `D.(Intercept)` = 11.1479626, `D.noise` = -0.5586893986, `mu` = 9.935279256),
    se.linked = c(`sigma_link` = 0.05708748926, `lambda0_link` = 0.2145093773, `sigma.toa_link` = 0.1019811062, `D.(Intercept)_link` = 4.747942207, `D.noise_link` = 0.449520798, `mu_link` = 0.1488572442),
    se.fitted = c(`sigma` = 0.1405154867, `lambda0` = 1.046461228, `sigma.toa` = 0.0002045132741, `D.(Intercept)` = 4.747942207, `D.noise` = 0.449520798, `mu` = 1.47893829)),
  ind_ss_log = list(
    linked = c(`b0.ss_link` = 1.67814029, `b1.ss_link` = -2.572589329, `sigma.ss_link` = 2.323011794, `D_link` = 4.517961021, `mu_link` = 2.161141087),
    fitted = c(`b0.ss` = 5.355586866, `b1.ss` = 0.07633762603, `sigma.ss` = 10.20636754, `D` = 91.64853787, `mu` = 8.681037831),
    se.linked = c(`b0.ss_link` = 0.001578161572, `b1.ss_link` = 0.01570423311, `sigma.ss_link` = 0.02537620096, `D_link` = 0.1740843391, `mu_link` = 0.06342224379),
    se.fitted = c(`b0.ss` = 0.008451981386, `b1.ss` = 0.001198823874, `sigma.ss` = 0.2589988337, `D` = 15.95457514, `mu` = 0.5505708977)),
  ind_ss_sp = list(
    linked = c(`b0.ss_link` = 4.625058463, `b1.ss_link` = -1.745818201, `sigma.ss_link` = 2.304772633, `D_link` = 4.235713518, `mu_link` = 2.247738099),
    fitted = c(`b0.ss` = 102.0087366, `b1.ss` = 0.1745021528, `sigma.ss` = 10.02189935, `D` = 69.11097306, `mu` = 9.466299774),
    se.linked = c(`b0.ss_link` = 0.01656569356, `b1.ss_link` = 0.550896673, `sigma.ss_link` = 0.02796464379, `D_link` = 0.2000028191, `mu_link` = 0.06543292427),
    se.fitted = c(`b0.ss` = 1.689845472, `b1.ss` = 0.09613265539, `sigma.ss` = 0.2802588453, `D` = 13.82238944, `mu` = 0.6194076763))
)

expect.fit.matches <- function(fit, expected, tolerance = 1e-4){
  expect_equal(coef(fit, types = 'linked'), expected$linked, tolerance = tolerance)
  expect_equal(coef(fit, types = 'fitted'), expected$fitted, tolerance = tolerance)
  expect_equal(stdEr(fit, types = 'linked'), expected$se.linked, tolerance = tolerance)
  expect_equal(stdEr(fit, types = 'fitted'), expected$se.fitted, tolerance = tolerance)
}

for(data.name in names(expected.fast)){
  test_that(paste0("fast fit on coarse mask -- ", data.name), {
    fit <- fit.demo.coarse(data.name)
    expect.fit.matches(fit, expected.fast[[data.name]])
  })
}

test_that("gradient-free fit agrees with the gradient-based one", {
  fit <- fit.demo.coarse('ihd_ext', gradient.free = TRUE)
  expect.fit.matches(fit, expected.fast[['ihd_ext']])
})

test_that("estimates, standard errors and intervals with new covariate values", {
  fit <- fit.demo.coarse('ihd_ext')
  new.data <- data.frame(brand = 'sony', noise = 7.7)
  pars <- c('g0', 'sigma', 'D')
  expect_equal(coef(fit, types = 'fitted', new.covariates = new.data)[pars], c(g0 = 1, sigma = 5.73517236, D = 831.7528031), tolerance = 1e-4)
  expect_equal(stdEr(fit, types = 'fitted', new.covariates = new.data)[c('sigma', 'D')], c(sigma = 0.9169049017, D = 399.8796185), tolerance = 1e-4)
  ci <- confint(fit, types = 'fitted', new.covariates = new.data)[pars, ]
  expect_equal(ci[, 1], c(g0 = 1.000000, sigma = 4.192388, D = 324.163666), tolerance = 1e-4)
  expect_equal(ci[, 2], c(g0 = 1.000000, sigma = 7.845696, D = 2134.146417), tolerance = 1e-4)
  expect_equal(colnames(ci), c('2.5 %', '97.5 %'))
  ## The detection function intercept g0 is fixed, so it has no standard error.
  expect_error(stdEr(fit, par = 'g0'))
})

test_that("confidence intervals are consistent with estimates and standard errors", {
  fit <- fit.demo.coarse('hhn_cue')
  ci <- confint(fit, types = 'linked', level = 0.9)
  est <- coef(fit, types = 'linked')[rownames(ci)]
  se <- stdEr(fit, types = 'linked')[rownames(ci)]
  expect_equal(unname(ci[, 1]), unname(est + stats::qnorm(0.05) * se), tolerance = 1e-6)
  expect_equal(unname(ci[, 2]), unname(est + stats::qnorm(0.95) * se), tolerance = 1e-6)
  expect_equal(colnames(ci), c('5 %', '95 %'))
})
