## Tests for p.dot(), which uses the detection probabilities reported by TMB at the masks used to
## fit the model, and otherwise calculates them in R.

test_that("detection probabilities from TMB match those calculated in R", {
  ## Session- and detector-level covariates over two sessions, detector- and mask-level
  ## covariates, signal strengths with an animal model, a cue-rate model, and a gradient-free fit.
  fits <- list(fit.demo.coarse('mul_ses_ext', nx = 15), fit.demo.coarse('ihd_ext'),
               fit.demo.coarse('ind_ss'), fit.demo.coarse('hhn_cue'),
               fit.demo.coarse('ihd_ext', gradient.free = TRUE))
  for(fit in fits){
    tmb <- p.dot(fit)
    r <- p.dot(fit, mask = get.mask(fit))
    expect_false(is.null(fit$output.tmb$p_dot))
    expect_equal(lengths(tmb), sapply(get.mask(fit), nrow))
    expect_equal(tmb, r, tolerance = 1e-8)
    expect_true(all(unlist(tmb) >= 0 & unlist(tmb) <= 1) && sd(unlist(tmb)) > 0)
  }
  ## Fits from before detection probabilities were reported use the calculation in R.
  fit$output.tmb$p_dot <- NULL
  expect_equal(p.dot(fit), r)
})

test_that("detection probabilities are the probability of detection by at least one detector", {
  ## A halfnormal detection function with g0 = 1, calculated directly.
  dat <- read.acre(captures = ihd$capt, traps = ihd$traps,
                   control.mask = list(buffer = ihd$control.mask$buffer, nx = 16))
  fit <- fit.quietly(dat, list(detfn = 'hn', optim.opts = list(fix = list(g0 = 1))))
  mask <- get.mask(fit)[[1]]
  traps <- as.matrix(ihd$traps)
  dists <- as.matrix(dist(rbind(mask, traps)))[1:nrow(mask), nrow(mask) + 1:nrow(traps)]
  p <- exp(-dists^2/(2*coef(fit, types = 'fitted')[['sigma']]^2))
  expect_equal(p.dot(fit)[[1]], unname(1 - apply(1 - p, 1, prod)), tolerance = 1e-8)
})

test_that("detection probabilities at new points use new covariate values", {
  ## sigma depends on a mask-level covariate.
  dat <- read.acre(captures = ihd$capt, traps = ihd$traps, loc.cov = ihd$loc.cov,
                   control.mask = list(buffer = ihd$control.mask$buffer, nx = 16))
  fit <- fit.quietly(dat, list(detfn = 'hn', model = list(sigma = ~ noise, D = ~ s(x, k = 5)),
                               optim.opts = list(fix = list(g0 = 1))))
  expect_true(abs(coef(fit)[['sigma.noise_link']]) > 1e-3)
  mask <- get.mask(fit)[[1]]
  noise <- get.par.extend.data(fit)$mask$noise
  rows <- c(10, 3, 50)
  ## At some of the fitted mask points, with their covariates, in another order.
  expect_equal(p.dot(fit, mask = mask[rows, ], newdata = data.frame(noise = noise[rows]))[[1]],
               p.dot(fit)[[1]][rows], tolerance = 1e-8)
  ## With other covariate values, sigma changes accordingly.
  new <- p.dot(fit, mask = mask[rows, ], newdata = data.frame(noise = noise[rows] + 1))[[1]]
  sigma <- exp(coef(fit)[['sigma.(Intercept)_link']] + coef(fit)[['sigma.noise_link']]*(noise[rows] + 1))
  traps <- as.matrix(ihd$traps)
  dists <- as.matrix(dist(rbind(mask[rows, ], traps)))[1:3, 3 + 1:nrow(traps)]
  expect_equal(new, unname(1 - apply(1 - exp(-dists^2/(2*sigma^2)), 1, prod)), tolerance = 1e-8)
  ## New points need their covariates.
  expect_error(p.dot(fit, mask = mask[rows, ]), "noise")
})
