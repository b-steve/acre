## Tests for the effects of terms in the model for D, plotted by plot(fit, type = "effect").

## The ihd demo data on a coarse mask, and a fit to it with a halfnormal detection function.
ihd.effect.fit <- function(model){
  dat <- read.acre(captures = ihd$capt, traps = ihd$traps, loc.cov = ihd$loc.cov,
                   control.mask = list(buffer = ihd$control.mask$buffer, nx = 16))
  fit.quietly(dat, list(model = list(D = model), detfn = 'hn', optim.opts = list(fix = list(g0 = 1))))
}

## Plots to a null device, returning what the plot returns.
quiet.plot <- function(...){
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  plot(...)
}

test_that("smooth effects are centred, match predictions, and have standard errors", {
  x <- seq(-20, 20, length.out = 15)
  for(model in list(~ y + s(x, k = 6), ~ y + s(x, k = 5, fx = TRUE))){
    fit <- ihd.effect.fit(model)
    eff <- term.effect(fit, "x", x)
    expect_equal(term.effect(fit, "s(x)", x), eff)
    ## The effect differs from the linear predictor only by a constant when y is held fixed.
    p <- predict(fit, newdata = data.frame(x = x, y = 3), type = "link", realnames = "D")$D$Estimate
    expect_equal(eff$est - p, rep(mean(eff$est - p), length(x)), tolerance = 1e-8)
    expect_true(sd(eff$est) > 0.01)
    ## Like mgcv's, the effect sums to zero over the data used to set up the model.
    mask <- get.mask(fit)[[1]]
    expect_equal(mean(term.effect(fit, "x", mask[, "x"])$est), 0, tolerance = 1e-8)
    ## Standard errors from the covariance of the smooth's coefficients only, following mgcv's
    ## recipe for smooth2random() predictions.
    s <- get.gam(fit, 'D')$gam_mask$smooth_info$smooths[[1]]
    if(is.null(s$re)){
      X <- mgcv::PredictMat(s$sm, data.frame(x = x))
    } else {
      r <- smooth2random.predict(s$sm, s$re, data.frame(x = x))
      X <- cbind(r$Xf, r$rand[[1]])
    }
    V <- vcov(fit, types = 'linked', pars = 'D', re = TRUE)
    keep <- grepl("^D\\.s\\(x\\)\\.", rownames(V))
    expect_equal(sum(keep), ncol(X))
    expect_equal(eff$se, sqrt(diag(X %*% V[keep, keep] %*% t(X))), tolerance = 1e-8)
    expect_true(all(eff$se > 0))
  }
})

test_that("linear effects are centred, and their standard errors are from their coefficients", {
  fit <- ihd.effect.fit(~ y + noise)
  mask_noise <- get.par.extend.data(fit)$mask$noise
  noise <- seq(5, 12, length.out = 10)
  eff <- term.effect(fit, "noise", noise)
  b <- coef(fit)[["D.noise_link"]]
  expect_equal(eff$est, b*(noise - mean(mask_noise)), tolerance = 1e-8)
  expect_equal(eff$se, abs(noise - mean(mask_noise))*stdEr(fit)[["D.noise_link"]], tolerance = 1e-8)
  ## A quadratic: both terms are included, with the covariance of their coefficients.
  fit <- ihd.effect.fit(~ y + noise + I(noise^2))
  eff <- term.effect(fit, "noise", noise)
  p <- predict(fit, newdata = data.frame(noise = noise, y = 3), type = "link", realnames = "D")$D$Estimate
  expect_equal(eff$est - p, rep(mean(eff$est - p), length(noise)), tolerance = 1e-8)
  expect_equal(mean(term.effect(fit, "noise", mask_noise)$est), 0, tolerance = 1e-8)
  X <- cbind(noise - mean(mask_noise), noise^2 - mean(mask_noise^2))
  pars <- c("D.noise_link", "D.I(noise^2)_link")
  V <- vcov(fit, types = 'linked', pars = 'D')[pars, pars]
  expect_equal(eff$se, sqrt(diag(X %*% V %*% t(X))), tolerance = 1e-8)
  expect_error(term.effect(fit, "forest_volumn", 1), "not a term or covariate")
  expect_error(term.effect(ihd.effect.fit(~ forest_volumn + noise:y), "noise", 1), "interactions")
})

test_that("effects of two covariates are surfaces", {
  fit <- ihd.effect.fit(~ noise + s(x, y, k = 8))
  mask <- get.mask(fit)[[1]]
  points <- data.frame(x = mask[c(1, 40, 90), "x"], y = mask[c(1, 40, 90), "y"])
  eff <- term.effect(fit, "s(x,y)", points)
  expect_equal(term.effect(fit, "x", points), eff)
  p <- predict(fit, newdata = cbind(points, noise = 8), type = "link", realnames = "D")$D$Estimate
  expect_equal(eff$est - p, rep(mean(eff$est - p), nrow(points)), tolerance = 1e-8)
  ## By default, on a grid of pixels within the mask's buffer of the detectors.
  out <- quiet.plot(fit, type = "effect", term = "s(x,y)", x.pixels = 10, y.pixels = 12, se = TRUE, ask = FALSE)
  grid <- out[["s(x,y)"]]
  expect_equal(grid$est, term.effect(fit, "s(x,y)", grid[, c("x", "y")])$est)
  traps <- as.matrix(ihd$traps)
  buffer <- attr(mask, "buffer")
  near <- function(p) apply(outer(p$x, traps[, 1], "-")^2 + outer(p$y, traps[, 2], "-")^2, 1, min) <= buffer^2
  all_pixels <- expand.grid(x = seq(min(traps[, 1]) - buffer, max(traps[, 1]) + buffer, length.out = 10),
                            y = seq(min(traps[, 2]) - buffer, max(traps[, 2]) + buffer, length.out = 12))
  expect_equal(grid[, c("x", "y")], all_pixels[near(all_pixels), ], ignore_attr = TRUE)
  expect_true(sum(near(all_pixels)) < 120)
  ## Including unsuitable habitat, which is not in the mask.
  dat <- read.acre(captures = ihd$capt, traps = ihd$traps, loc.cov = ihd$loc.cov,
                   control.mask = list(buffer = ihd$control.mask$buffer, nx = 16),
                   unsuitable = forest_volumn == "low")
  fit_u <- fit.quietly(dat, list(model = list(D = ~ s(x, y, k = 8)), detfn = 'hn',
                                 optim.opts = list(fix = list(g0 = 1))))
  grid_u <- quiet.plot(fit_u, type = "effect", term = "s(x,y)", x.pixels = 10, y.pixels = 12)[[1]]
  expect_equal(grid_u[, c("x", "y")], grid[, c("x", "y")])
  low <- location.cov.to.mask(list(as.matrix(grid_u[, c("x", "y")])), ihd$loc.cov)$forest_volumn == "low"
  expect_true(any(low))
  ## At points provided.
  out <- quiet.plot(fit, type = "effect", term = "s(x,y)", newdata = points)
  expect_equal(out[["s(x,y)"]], eff)
})

test_that("effect plots show what they return, over the mask points or detectable ones", {
  fit <- ihd.effect.fit(~ y + s(x, k = 6))
  mask <- get.mask(fit)[[1]]
  out <- quiet.plot(fit, type = "effect")
  expect_equal(names(out), "s(x)")
  expect_equal(out[[1]]$est, term.effect(fit, "x", out[[1]]$x)$est)
  expect_equal(range(out[[1]]$x), range(mask[, "x"]))
  ## Only the detectable mask points.
  detectable <- p.dot(fit)[[1]] > 0.1
  expect_true(any(detectable) && !all(detectable))
  out <- quiet.plot(fit, type = "effect", p.dot.min = 0.1)
  expect_equal(range(out[[1]]$x), range(mask[detectable, "x"]))
  out <- quiet.plot(fit, type = "effect", xlim = c(-5, 5))
  expect_equal(range(out[[1]]$x), c(-5, 5))
  ## Adding to an existing plot, e.g., a linear effect and a smooth one.
  fit_lin <- ihd.effect.fit(~ y + x)
  grDevices::pdf(NULL)
  plot(fit_lin, type = "effect", term = "x")
  out <- plot(fit, type = "effect", term = "x", add = TRUE, col = "red")
  grDevices::dev.off()
  expect_equal(out[[1]]$est, term.effect(fit, "x", out[[1]]$x)$est)
  expect_error(quiet.plot(ihd.effect.fit(~ y), type = "effect"), "no smooths")
})
