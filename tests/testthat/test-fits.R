
test_that("joint bearing/dist fitting -- half normal", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[1])$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(36.04365, 1.607376, 3.876604, 1.892325, 7.835287)

  pars_link_names = c("g0_link", "sigma_link", "kappa_link", "alpha_link", "D_link")
  pars_names = c('g0', 'sigma', 'kappa', 'alpha', 'D')
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)

  #test fitted estimations
  pars_est_values = c(1.000000, 4.989701, 48.260046, 6.634774, 2528.260710)
  relative.error = max(abs((coef(fit, types = 'fitted')[pars_names] - pars_est_values)/pars_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #g0 is fixed, so it has no std error estimation
  pars_link_std_values = c(0.03049075, 0.1397999, 0.1059173, 0.101818)
  
  #test linked std errors
  pars_link_names_unfixed = pars_link_names[-1]
  pars_names_unfixed = pars_names[-1]
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names_unfixed] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #expect an error here as g0 is fixed, no estimated std error
  expect_error(stdEr(fit, par = 'g0'))
  
  #test fitted std errors
  pars_std_values = c(0.1521397, 6.7467492, 0.7027375, 257.4224243)
  relative.error = max(abs((stdEr(fit, types = 'fitted')[pars_names_unfixed] - pars_std_values)/pars_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(36.043653, 1.547615, 3.602601, 1.684730, 7.635727, 36.043653, 1.667137,
                     4.150607, 2.099919, 8.034846), ncol = 2)
  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  #test another confidence level
  conf_90 = matrix(c(36.043653, 1.557223, 3.646654, 1.718106, 7.667811, 36.043653, 1.657529,
                     4.106554, 2.066543,  8.002763), ncol = 2)

  
  o = confint(fit, types = 'linked', level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  #test fitted confident interval
  conf_95_fitted = matrix(c(1.000000, 4.700248, 36.693560, 5.390997, 2070.876621, 1.000000, 5.296980,
                            63.472502, 8.165506, 3086.664919), ncol = 2)
  o = confint(fit, types = 'fitted', level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted)/conf_95_fitted))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  ############################################################################################################
  #since this model has no parameter been extended, skip the test for 
  
})


test_that("heterogeneous density model with sigma extended -- half normal", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[5])$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(36.0436534, 2.0142606, -0.2688631, 5.7910711, 0.1208786)
  
  pars_link_names = c("g0_link", "sigma.(Intercept)_link", "sigma.brandsony_link", "D.(Intercept)_link", "D.noise_link")
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #g0 is fixed, so it has no std error estimation
  pars_link_std_values = c(0.13479140, 0.06918262, 2.06221583, 0.21047749)
  
  #test linked std errors
  pars_link_names_unfixed = pars_link_names[-1]
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names_unfixed] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #expect an error here as g0 is fixed, no estimated std error
  expect_error(stdEr(fit, par = 'g0'))
  
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(36.0436534, 1.7500743, -0.4044586, 1.7492023, -0.2916497, 36.0436534, 2.2784469,
                     -0.1332677, 9.8329399, 0.5334069), ncol = 2)

  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  
  #test another confidence level
  conf_90 = matrix(c(36.0436534, 1.7925485, -0.3826584, 2.3990279, -0.2253260, 36.0436534, 2.2359727,
                     -0.1550679, 9.1831143,  0.4670833), ncol = 2)

  
  o = confint(fit, types = 'linked', level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  

  
  ############################################################################################################
  #test for new_data provided
  new_data = data.frame(brand = 'sony', noise = 7.7)
  pars_names_og = c('g0', 'sigma', 'D')
  pars_names_og_unfixed = pars_names_og[-1]
  
  expected_values = c(1, 5.728178, 830.341102)
  o = coef(fit, types = 'fitted', new.covariates = new_data)[pars_names_og]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  expected_values = c(0.9222248, 400.2406206)
  o = stdEr(fit, types = 'fitted', new.covariates = new_data)[pars_names_og_unfixed]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  expected_values = matrix(c(1, 4.178061, 322.820254, 1, 7.853408, 2135.759253),
                           ncol = 2)
  o = confint(fit, types = 'fitted', new.covariates = new_data)[pars_names_og,]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
})



test_that("joint bearing/dist model with 2 sessions and g0 & sigma extended -- half normal", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[7])$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(1.5378647, -1.0670451, 1.0895859, 0.1642574, 3.5215606, 1.0289551, 7.8714508)
  
  pars_link_names = c("g0.(Intercept)_link", "g0.weathersunny_link", "sigma.(Intercept)_link", "sigma.brandsony_link",
                      "kappa_link", "alpha_link", "D_link")
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  #test linked std errors
  pars_link_std_values = c(1.58088441, 1.36630628, 0.07671678, 0.07736679, 0.51800547, 0.21100953, 0.17920083)
 
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(-1.56061182, -3.74495623, 0.93922376, 0.01262129, 2.50628848, 0.61538402, 7.52022367,
                     4.6363412, 1.6108660, 1.2399480, 0.3158935, 4.5368326, 1.4425262, 8.2226780), ncol = 2)
  
  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  
  #test another confidence level
  conf_90 = matrix(c(-1.06245877, -3.31441897, 0.96339802, 0.03700037, 2.66951737, 0.68187531, 7.57669171,
                     4.1381881, 1.1803287, 1.2157738, 0.2915144, 4.3736037, 1.3760349, 8.1662100), ncol = 2)

  o = confint(fit, types = 'linked', level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  
  #test fitted confident interval
  pars_names = c("g0.(Intercept)", "g0.weathersunny", "sigma.(Intercept)", "sigma.brandsony",
                 "kappa", "alpha", "D")
  
  conf_95_fitted = matrix(c(-1.56061182, -3.74495623, 0.93922376, 0.01262129, 12.25934474, 1.85036704, 1844.97990856,
                            4.6363412, 1.6108660, 1.2399480, 0.3158935, 93.3945155, 4.2313715, 3724.4632673), ncol = 2)
  o = confint(fit, types = 'fitted', level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted)/conf_95_fitted))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  #test for new_data provided
  new_data = data.frame(brand = 'sony', weather = 'sunny')
  pars_names_og = c('g0', 'sigma', 'kappa', 'alpha', 'D')
  
  expected_values = c(0.6155777, 3.5037832, 33.8371920, 2.7981405, 2621.3660367)
  o = coef(fit, types = 'fitted', new.covariates = new_data)[pars_names_og]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  expected_values = c(0.1210091, 0.2763892, 17.5278506, 0.5904343, 469.7509821)
  o = stdEr(fit, types = 'fitted', new.covariates = new_data)[pars_names_og]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  expected_values = matrix(c(0.370184, 3.001870, 12.259345, 1.850367, 1844.979909,
                             0.8135231, 4.0896168, 93.3945155, 4.2313715, 3724.4632673),
                           ncol = 2)
  o = confint(fit, types = 'fitted', new.covariates = new_data)[pars_names_og,]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
})


test_that("cue rate included -- hazard half normal", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[9])$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(1.405026, 1.519666, 6.203402)
  
  pars_link_names = c("sigma_link", "lambda0_link", "D_link")
  pars_names = c('sigma', 'lambda0', 'D')
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #test fitted estimations
  pars_est_values = c(4.075632, 4.570699, 494.428405)
  relative.error = max(abs((coef(fit, types = 'fitted')[pars_names] - pars_est_values)/pars_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  #test linked std errors
  pars_link_std_values = c(0.07907883, 0.19694115, 0.13699999)
  
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  #test fitted std errors
  pars_std_values = c(0.3222962, 0.9001587, 67.7366874)
  relative.error = max(abs((stdEr(fit, types = 'fitted')[pars_names] - pars_std_values)/pars_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(1.250034, 1.133669, 5.934887, 1.560017, 1.905664, 6.471917), ncol = 2)
  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  #test another confidence level
  conf_90 = matrix(c(1.274953, 1.195727, 5.978057, 1.535099, 1.843605, 6.428747), ncol = 2)
  
  
  o = confint(fit, types = 'linked', level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  #test fitted confident interval
  conf_95_fitted = matrix(c(3.490462, 3.107034, 377.997397, 4.758904, 6.723868, 646.722570), ncol = 2)
  o = confint(fit, types = 'fitted', level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted)/conf_95_fitted))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  ############################################################################################################
  #since this model has no parameter been extended, skip the test for 
  
})


test_that("Signal strength & toa model", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[12])$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(4.497227, 1.357671, 2.227983, -6.400474, 7.814611)
  
  pars_link_names = c("b0.ss_link", "b1.ss_link", "sigma.ss_link", "sigma.toa_link", "D_link")
  pars_names = c('b0.ss', 'b1.ss', 'sigma.ss', 'sigma.toa', 'D')
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #test fitted estimations
  pars_est_values = c(8.976785e+01, 3.887131e+00, 9.281129e+00, 1.660770e-03, 2.476523e+03)
  relative.error = max(abs((coef(fit, types = 'fitted')[pars_names] - pars_est_values)/pars_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test linked std errors
  pars_link_std_values = c(0.01398521, 0.05075748, 0.04960075, 0.08670580, 0.09944047)
  
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  #test fitted std errors
  pars_std_values = c(1.255422e+00, 1.973010e-01, 4.603510e-01, 1.439984e-04, 2.462666e+02)
  relative.error = max(abs((stdEr(fit, types = 'fitted')[pars_names] - pars_std_values)/pars_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(4.469816, 1.258189, 2.130767, -6.570414, 7.619711,
                     4.524637, 1.457154, 2.325199, -6.230534, 8.009511), ncol = 2)
  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  #test another confidence level
  conf_90 = matrix(c(4.474223, 1.274183, 2.146397, -6.543093, 7.651046,
                     4.520231, 1.44116, 2.309569, -6.257856, 7.978176), ncol = 2)
  
  
  o = confint(fit, types = 'linked', level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  #test fitted confident interval
  conf_95_fitted = matrix(c(87.34069, 3.519041, 8.421327, 0.001401217, 2037.974,
                            92.26247, 4.293723, 10.22871, 0.001968401, 3009.444), ncol = 2)
  o = confint(fit, types = 'fitted', level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted)/conf_95_fitted))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  ############################################################################################################
  #since this model has no parameter been extended, skip the test for 
  
})


test_that("heterogeneous density & toa model with individual identity -- hazard half normal", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[14])$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(0.7707158, 1.835558, -6.63785, 10.60063, -0.4870564, 2.165755)
  
  pars_link_names = c("sigma_link", "lambda0_link", "sigma.toa_link", "D.(Intercept)_link", "D.noise_link", "mu_link")
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test std error
  pars_link_std_values = c(0.07506491, 0.2700223, 0.1091684, 5.117405, 0.4782015, 0.1426855)
  
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(0.6235912, 1.306324, -6.851816, 0.5707017, -1.424314, 1.886096,
                     0.9178403, 2.364792, -6.423884, 20.63056, 0.4502012, 2.445413), ncol = 2)
  
  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  
  #test another confidence level
  conf_90 = matrix(c(0.647245, 1.391411, -6.817416, 2.183249, -1.273628, 1.931058,
                     0.8941865, 2.279705, -6.458284, 19.01801, 0.299515, 2.400452), ncol = 2)
  
  
  o = confint(fit, types = 'linked', level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  
  
  ############################################################################################################
  #test for new_data provided
  new_data = data.frame(noise = 7.7)
  pars_names_og = c('sigma', 'lambda0', 'sigma.toa', 'D', 'mu')
  
  expected_values = c(2.161313, 6.26863, 0.00130984, 944.1608, 8.721183)
  o = coef(fit, types = 'fitted', new.covariates = new_data)[pars_names_og]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  expected_values = c(0.1622387, 1.69267, 0.0001429931, 1372.344, 1.244387)
  o = stdEr(fit, types = 'fitted', new.covariates = new_data)[pars_names_og]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  expected_values = matrix(c(1.865616, 3.692574, 0.001057533, 54.67874, 6.593579,
                             2.503877, 10.64182, 0.001622343, 16303.22, 11.53532),
                           ncol = 2)
  o = confint(fit, types = 'fitted', new.covariates = new_data)[pars_names_og,]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
})


test_that("signal strength model with log link and individual identity", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[16])$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(1.676551, -2.532202, 2.287277, 4.529954, 2.162836)
  
  pars_link_names = c("b0.ss_link", "b1.ss_link", "sigma.ss_link", "D_link", "mu_link")
  pars_names = c('b0.ss', 'b1.ss', 'sigma.ss', 'D', 'mu')
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #test fitted estimations
  pars_est_values = c(5.3470820, 0.0794838, 9.8480851, 92.7542590, 8.6957638)
  relative.error = max(abs((coef(fit, types = 'fitted')[pars_names] - pars_est_values)/pars_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #g0 is fixed, so it has no std error estimation
  pars_link_std_values = c(0.001742142, 0.015412422, 0.025178517, 0.174155514, 0.063878707)
  
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  #test fitted std errors
  pars_std_values = c(0.009315378, 0.001225038, 0.247960173, 16.153665696, 0.555474145)
  relative.error = max(abs((stdEr(fit, types = 'fitted')[pars_names] - pars_std_values)/pars_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(1.673136, -2.562410, 2.237928, 4.188615, 2.037636,
                     1.679966, -2.501994, 2.336626, 4.871292, 2.288036), ncol = 2)
  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  #test another confidence level
  conf_90 = matrix(c(1.673685, -2.557553, 2.245862, 4.243493, 2.057765,
                     1.679417, -2.506851, 2.328692, 4.816414, 2.267907), ncol = 2)
  
  
  o = confint(fit, types = 'linked', level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  #test fitted confident interval
  conf_95_fitted = matrix(c(5.32885530, 0.07711868, 9.37388885, 65.93141820, 7.67245027,
                            5.36537098, 0.08192147, 10.34626944, 130.48942065, 9.85556188), ncol = 2)
  
  o = confint(fit, types = 'fitted', level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted)/conf_95_fitted))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  ############################################################################################################
  #since this model has no parameter been extended, skip the test for 
  
})


test_that("signal strength model with spherical link and individual identity", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[17])$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(4.627961, -1.679443, 2.305293, 4.234750, 2.247667)
  
  pars_link_names = c("b0.ss_link", "b1.ss_link", "sigma.ss_link", "D_link", "mu_link")
  pars_names = c('b0.ss', 'b1.ss', 'sigma.ss', 'D', 'mu')
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #test fitted estimations
  pars_est_values = c(102.3052278, 0.1864778, 10.0271191, 69.0444275, 9.4656292)
  relative.error = max(abs((coef(fit, types = 'fitted')[pars_names] - pars_est_values)/pars_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #g0 is fixed, so it has no std error estimation
  pars_link_std_values = c(0.02051216, 0.61022968, 0.02797164, 0.19999997, 0.06543341)
  
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  #test fitted std errors
  pars_std_values = c(2.0985008, 0.1137943, 0.2804750, 13.8088834, 0.6193684)
  relative.error = max(abs((stdEr(fit, types = 'fitted')[pars_names] - pars_std_values)/pars_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(4.58776, -2.875357, 2.250471, 3.84273, 2.119423,
                     4.668163, -0.4834706, 2.360118, 4.626725, 2.375917), ncol = 2)
  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  #test another confidence level
  conf_90 = matrix(c(4.594221, -2.683182, 2.259284, 3.905779, 2.140039,
                     4.6617003, -0.6757045, 2.3513026, 4.5637209, 2.3552956), ncol = 2)
  
  
  o = confint(fit, types = 'linked', level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  #test fitted confident interval
  conf_95_fitted = matrix(c(98.27403, 0.05639598, 9.492202, 46.65267, 8.326332,
                            106.5019, 0.6166396, 10.5922, 102.1788, 10.76088), ncol = 2)
  
  o = confint(fit, types = 'fitted', level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted)/conf_95_fitted))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  ############################################################################################################
  #since this model has no parameter been extended, skip the test for 
  
})


test_that("heterogeneous density model with sigma extended -- half normal - no gradient", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[5], gradient_free = TRUE)$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(36.0436534, 2.0142606, -0.2688632, 5.7910710, 0.1208786)
  
  pars_link_names = c("g0_link", "sigma.(Intercept)_link", "sigma.brandsony_link", "D.(Intercept)_link", "D.noise_link")
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  pars_link_std_values = c(0.13479106, 0.06918256, 2.06221107, 0.21047694)
  
  #test linked std errors
  pars_link_names_unfixed = pars_link_names[-1]
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names_unfixed] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #expect an error here as g0 is fixed, no estimated std error
  expect_error(stdEr(fit, par = 'g0'))
  
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(36.0436534, 1.7500750, -0.4044585, 1.7492116, -0.2916486, 36.0436534, 2.2784462,
                     -0.1332678, 9.8329304, 0.5334059), ncol = 2)
  
  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  
  #test another confidence level
  conf_90 = matrix(c(36.0436534, 1.7925490, -0.3826583, 2.3990357, -0.2253251, 36.0436534, 2.2359722,
                     -0.1550680, 9.1831064,  0.4670824), ncol = 2)
  
  
  o = confint(fit, types = 'linked', level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  
  
  ############################################################################################################
  #test for new_data provided
  new_data = data.frame(brand = 'sony', noise = 7.7)
  pars_names_og = c('g0', 'sigma', 'D')
  pars_names_og_unfixed = pars_names_og[-1]
  
  expected_values = c(1, 5.728178, 830.341085)
  o = coef(fit, types = 'fitted', new.covariates = new_data)[pars_names_og]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  expected_values = c(0.9222219, 400.2402285)
  o = stdEr(fit, types = 'fitted', new.covariates = new_data)[pars_names_og_unfixed]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)

  expected_values = matrix(c(1, 4.178065, 322.820540, 1, 7.8534, 2135.7573),
                           ncol = 2)
  o = confint(fit, types = 'fitted', new.covariates = new_data)[pars_names_og,]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
})


test_that("Signal strength & toa model - no gradient", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[12], gradient_free = TRUE)$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(4.497227, 1.357671, 2.227983, -6.400474, 7.814611)
  
  pars_link_names = c("b0.ss_link", "b1.ss_link", "sigma.ss_link", "sigma.toa_link", "D_link")
  pars_names = c('b0.ss', 'b1.ss', 'sigma.ss', 'sigma.toa', 'D')
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #test fitted estimations
  pars_est_values = c(89.76786, 3.887131, 9.281129, 0.001660769, 2476.524)
  relative.error = max(abs((coef(fit, types = 'fitted')[pars_names] - pars_est_values)/pars_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test linked std errors
  pars_link_std_values = c(0.01398517, 0.05075751, 0.04960074, 0.08670579, 0.09944046)
  
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  #test fitted std errors
  pars_std_values = c(1.255419, 0.1973011, 0.4603509, 0.0001439983, 246.2667)
  relative.error = max(abs((stdEr(fit, types = 'fitted')[pars_names] - pars_std_values)/pars_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(4.469817, 1.258188, 2.130768, -6.570415, 7.619711,
                     4.524637, 1.457154, 2.325199, -6.230534, 8.009511), ncol = 2)
  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  #test another confidence level
  conf_90 = matrix(c(4.474223, 1.274183, 2.146397, -6.543093, 7.651046,
                     4.520231, 1.44116, 2.309569, -6.257856, 7.978176), ncol = 2)
  
  
  o = confint(fit, types = 'linked', level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  #test fitted confident interval
  conf_95_fitted = matrix(c(87.3407, 3.519041, 8.421328, 0.001401216, 2037.974,
                            92.26246, 4.293723, 10.22871, 0.0019684, 3009.445), ncol = 2)
  o = confint(fit, types = 'fitted', level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted)/conf_95_fitted))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  ############################################################################################################
  #since this model has no parameter been extended, skip the test for 
  
})


test_that("heterogeneous density & toa model with individual identity -- hazard half normal - no gradient", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[14], gradient_free = TRUE)$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(0.7707153, 1.835559, -6.637848, 10.60045, -0.4870398, 2.165754)
  
  pars_link_names = c("sigma_link", "lambda0_link", "sigma.toa_link", "D.(Intercept)_link", "D.noise_link", "mu_link")
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test std error
  pars_link_std_values = c(0.07506419, 0.2700215, 0.1091684, 5.117405, 0.4782013, 0.1426855)
  
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(0.6235922, 1.306327, -6.851815, 0.5705239, -1.424297, 1.886096,
                     0.9178384, 2.364792, -6.423882, 20.63038, 0.4502174, 2.445412), ncol = 2)
  
  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  
  #test another confidence level
  conf_90 = matrix(c(0.6472457, 1.391414, -6.817414, 2.183071, -1.273611, 1.931057,
                     0.8941849, 2.279705, -6.458282, 19.01784, 0.2995312, 2.400451), ncol = 2)
  
  
  o = confint(fit, types = 'linked', level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  
  
  ############################################################################################################
  #test for new_data provided
  new_data = data.frame(noise = 7.7)
  pars_names_og = c('sigma', 'lambda0', 'sigma.toa', 'D', 'mu')
  
  expected_values = c(2.161312, 6.26864, 0.001309842, 944.1142, 8.721175)
  o = coef(fit, types = 'fitted', new.covariates = new_data)[pars_names_og]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  expected_values = c(0.1622371, 1.692667, 0.0001429934, 1372.279, 1.244385)
  o = stdEr(fit, types = 'fitted', new.covariates = new_data)[pars_names_og]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  expected_values = matrix(c(1.865618, 3.692587, 0.001057535, 54.67584, 6.593574,
                             2.503872, 10.64182, 0.001622346, 16302.48, 11.53531),
                           ncol = 2)
  o = confint(fit, types = 'fitted', new.covariates = new_data)[pars_names_og,]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
})


test_that("signal strength model with spherical link and individual identity - no gradient", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[17], gradient_free = TRUE)$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(4.627962, -1.679417, 2.305293, 4.234735, 2.247669)
  
  pars_link_names = c("b0.ss_link", "b1.ss_link", "sigma.ss_link", "D_link", "mu_link")
  pars_names = c('b0.ss', 'b1.ss', 'sigma.ss', 'D', 'mu')
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #test fitted estimations
  pars_est_values = c(102.3053100, 0.1864827, 10.0271110, 69.0433944, 9.4656476)
  relative.error = max(abs((coef(fit, types = 'fitted')[pars_names] - pars_est_values)/pars_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  pars_link_std_values = c(0.02050358, 0.61003840, 0.02797149, 0.20000145, 0.06543333)
  
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  #test fitted std errors
  pars_std_values = c(2.0976247, 0.1137616, 0.2804732, 13.8087791, 0.6193689)
  relative.error = max(abs((stdEr(fit, types = 'fitted')[pars_names] - pars_std_values)/pars_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(4.587775, -2.875070, 2.250469, 3.842740, 2.119422,
                     4.6681478, -0.4837633, 2.3601156, 4.6267309, 2.3759162), ncol = 2)
  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  #test another confidence level
  conf_90 = matrix(c(4.594236, -2.682840, 2.259284, 3.905762, 2.140041,
                     4.6616870, -0.6759927, 2.3513015, 4.5637083, 2.3552975), ncol = 2)
  
  
  o = confint(fit, types = 'linked', level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  #test fitted confident interval
  conf_95_fitted = matrix(c(98.2754, 0.05640565, 9.492199, 46.6532, 8.326308,
                            106.5001, 0.6164695, 10.59219, 102.1796, 10.76085), ncol = 2)
  
  o = confint(fit, types = 'fitted', level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted)/conf_95_fitted))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  ############################################################################################################
  #since this model has no parameter been extended, skip the test for 
  
})

