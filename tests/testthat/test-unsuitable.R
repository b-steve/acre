## Tests for unsuitable habitat, which is removed from the mask.

## The ihd demo data on a coarse mask, with unsuitable habitat removed by read.acre().
ihd.unsuitable <- function(...){
  read.acre(captures = ihd$capt, traps = ihd$traps, control.mask = list(buffer = ihd$control.mask$buffer, nx = 16),
            loc.cov = ihd$loc.cov, ...)
}

test_that("read.acre() removes unsuitable habitat from the mask and its covariates", {
  full <- ihd.unsuitable()
  dat <- ihd.unsuitable(unsuitable = (forest_volumn == "low") | (x > 20 & noise > 9))
  cov_full <- full$par.extend$data$mask
  bad <- (cov_full$forest_volumn == "low") | (full$mask[[1]][, "x"] > 20 & cov_full$noise > 9)
  expect_true(any(bad) && !all(bad))
  expect_equal(dat$mask[[1]], full$mask[[1]][!bad, ], ignore_attr = TRUE)
  expect_equal(attr(dat$mask[[1]], "area"), attr(full$mask[[1]], "area"))
  expect_equal(attr(dat$mask[[1]], "buffer"), attr(full$mask[[1]], "buffer"))
  cov <- dat$par.extend$data$mask
  expect_equal(cov$mask, seq_len(sum(!bad)))
  expect_equal(cov[, c("noise", "forest_volumn")], cov_full[!bad, c("noise", "forest_volumn")], ignore_attr = TRUE)
  ## Variables from the calling environment, and one-sided formulas, also work.
  cutoff <- 9
  dat_cutoff <- read.acre(captures = ihd$capt, traps = ihd$traps, loc.cov = ihd$loc.cov,
                          control.mask = list(buffer = ihd$control.mask$buffer, nx = 16), unsuitable = noise > cutoff)
  expect_equal(dat_cutoff$mask[[1]], full$mask[[1]][cov_full$noise <= 9, ], ignore_attr = TRUE)
  expect_equal(ihd.unsuitable(unsuitable = ~ noise > 9)$mask, dat_cutoff$mask)
})

test_that("fit.acre() removes unsuitable habitat, including read.acre()'s from a new mask", {
  model <- list(D = ~ noise)
  fix <- list(optim.opts = list(fix = list(g0 = 1)))
  ## The fit with a mask that has had unsuitable habitat removed by hand.
  full <- ihd.unsuitable()
  bad <- full$par.extend$data$mask$forest_volumn == "low"
  mask_hand <- full$mask[[1]][!bad, ]
  attr(mask_hand, "area") <- attr(full$mask[[1]], "area")
  attr(mask_hand, "buffer") <- attr(full$mask[[1]], "buffer")
  dat_hand <- read.acre(captures = ihd$capt, traps = ihd$traps, mask = mask_hand, loc.cov = ihd$loc.cov)
  fit_hand <- fit.quietly(dat_hand, c(list(model = model), fix))
  ## Removing it with read.acre() or fit.acre() gives the same fit.
  fit_read <- fit.quietly(ihd.unsuitable(unsuitable = forest_volumn == "low"), c(list(model = model), fix))
  fit_fit <- fit.quietly(full, c(list(model = model, unsuitable = quote(forest_volumn == "low")), fix))
  for(fit in list(fit_read, fit_fit)){
    expect_equal(get.mask(fit), get.mask(fit_hand))
    expect_equal(fit$loglik, fit_hand$loglik)
    expect_equal(coef(fit), coef(fit_hand))
  }
  expect_equal(length(fit_fit$unsuitable), 1)
  ## A new mask from fit.acre() has unsuitable habitat from both read.acre() and fit.acre() removed.
  control_new <- list(buffer = ihd$control.mask$buffer, nx = 12)
  fit_new <- fit.quietly(ihd.unsuitable(unsuitable = forest_volumn == "low"),
                         c(list(model = model, unsuitable = quote(noise > 10), control.mask = control_new), fix))
  mask_new <- create.mask(ihd$traps, buffer = control_new$buffer, nx = control_new$nx)
  cov_new <- location.cov.to.mask(list(mask_new), ihd$loc.cov)
  bad_new <- cov_new$forest_volumn == "low" | cov_new$noise > 10
  expect_true(any(bad_new) && !all(bad_new))
  expect_equal(get.mask(fit_new)[[1]], mask_new[!bad_new, ], ignore_attr = TRUE)
  expect_equal(get.par.extend.data(fit_new)$mask$noise, cov_new$noise[!bad_new])
  expect_equal(length(fit_new$unsuitable), 2)
  ## Predicted density away from the mask is zero in unsuitable habitat, and otherwise as if
  ## there were no unsuitable habitat.
  grid <- list(xlim = c(-30, 35), ylim = c(-30, 40), x.pixels = 20, y.pixels = 20, se.fit = TRUE)
  pred <- do.call(D.predict.for.plot, c(list(fit_new), grid))
  fit_all <- fit_new
  fit_all$unsuitable <- list()
  pred_all <- do.call(D.predict.for.plot, c(list(fit_all), grid))
  cov_grid <- location.cov.to.mask(list(as.matrix(pred[, c("x", "y")])), ihd$loc.cov)
  bad_grid <- cov_grid$forest_volumn == "low" | cov_grid$noise > 10
  expect_true(any(bad_grid) && !all(bad_grid))
  expect_equal(pred$est[bad_grid], rep(0, sum(bad_grid)))
  expect_equal(pred$std[bad_grid], rep(0, sum(bad_grid)))
  expect_equal(pred$est[!bad_grid], pred_all$est[!bad_grid])
  expect_true(all(pred_all$est > 0))
})

test_that("unsuitable habitat is checked", {
  expect_error(ihd.unsuitable(unsuitable = noise + 1), "must give TRUE or FALSE")
  expect_error(ihd.unsuitable(unsuitable = elevation > 1), "Cannot evaluate 'unsuitable'")
  expect_error(ihd.unsuitable(unsuitable = noise > 0), "All mask points for session 1")
  expect_error(ihd.unsuitable(unsuitable = noise ~ 1), "one-sided")
})
