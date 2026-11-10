## Tests for covariates attached to masks, as columns other than x and y.

## A coarse mask for the ihd demo data, with covariate values that interpolation would not give.
ihd.mask.cov <- function(){
  mask <- create.mask(ihd$traps, buffer = ihd$control.mask$buffer, nx = 16)
  df <- data.frame(x = mask[, "x"], y = mask[, "y"])
  df$noise <- round(5 + df$x/10 + (df$y/20)^2, 3)
  df$forest_volumn <- ifelse(df$x > df$y, "low", "high")
  attr(df, "area") <- attr(mask, "area")
  attr(df, "buffer") <- attr(mask, "buffer")
  df
}

test_that("covariates attached to a mask are used as they are", {
  mask <- ihd.mask.cov()
  dat <- read.acre(captures = ihd$capt, traps = ihd$traps, mask = mask)
  expect_equal(unname(dat$mask[[1]]), unname(as.matrix(mask[, c("x", "y")])), ignore_attr = TRUE)
  expect_equal(colnames(dat$mask[[1]]), c("x", "y"))
  expect_equal(attr(dat$mask[[1]], "area"), attr(mask, "area"))
  expect_equal(attr(dat$mask[[1]], "buffer"), attr(mask, "buffer"))
  cov <- dat$par.extend$data$mask
  expect_equal(cov$mask, seq_len(nrow(mask)))
  expect_equal(cov$noise, mask$noise)
  expect_equal(cov$forest_volumn, mask$forest_volumn)
  ## Factors are treated like characters.
  mask_f <- mask
  mask_f$forest_volumn <- factor(mask_f$forest_volumn)
  expect_equal(read.acre(captures = ihd$capt, traps = ihd$traps, mask = mask_f)$par.extend, dat$par.extend)
  ## The fit is the same as with the covariates provided at exactly the mask points in loc.cov.
  model <- list(D = ~ noise + forest_volumn)
  fix <- list(fix = list(g0 = 1))
  fit_att <- fit.quietly(dat, list(model = model, optim.opts = fix))
  dat_loc <- read.acre(captures = ihd$capt, traps = ihd$traps, mask = as.matrix(mask[, c("x", "y")]),
                       loc.cov = mask)
  expect_equal(dat_loc$par.extend$data$mask[, c("noise", "forest_volumn")], cov[, c("noise", "forest_volumn")])
  fit_loc <- fit.quietly(dat_loc, list(model = model, optim.opts = fix))
  expect_equal(fit_att$loglik, fit_loc$loglik)
  expect_equal(coef(fit_att), coef(fit_loc))
  ## Unsuitable habitat is removed with its covariates.
  dat_u <- read.acre(captures = ihd$capt, traps = ihd$traps, mask = mask, unsuitable = noise > 8)
  keep <- mask$noise <= 8
  expect_true(any(!keep) && any(keep))
  expect_equal(unname(dat_u$mask[[1]]), unname(as.matrix(mask[keep, c("x", "y")])), ignore_attr = TRUE)
  expect_equal(dat_u$par.extend$data$mask$noise, mask$noise[keep])
})

test_that("covariates attached to a mask are used in favour of those provided elsewhere", {
  mask <- ihd.mask.cov()
  loc.cov <- ihd$loc.cov
  loc.cov$elevation <- loc.cov$x
  ## From loc.cov: noise is ignored, but the other covariates are still interpolated.
  expect_warning(dat <- read.acre(captures = ihd$capt, traps = ihd$traps, mask = mask[, c("x", "y", "noise")],
                                  loc.cov = loc.cov),
                 "Covariate 'noise' in 'loc.cov' is ignored in favour of the covariate attached to the mask.")
  interp <- location.cov.to.mask(list(as.matrix(mask[, c("x", "y")])), loc.cov)
  cov <- dat$par.extend$data$mask
  expect_equal(cov$noise, mask$noise)
  expect_equal(cov$forest_volumn, interp$forest_volumn)
  expect_equal(cov$elevation, interp$elevation)
  expect_false("noise" %in% colnames(dat$arg.input$loc.cov))
  ## Several at once, and from dist.cov, session.cov and trap.cov.
  dist.cov <- list(noise = data.frame(x = 0, y = 0), road = data.frame(x = c(0, 10), y = c(0, 0)))
  expect_warning(expect_warning(expect_warning(
    dat <- read.acre(captures = ihd$capt, traps = ihd$traps, mask = mask, loc.cov = loc.cov, dist.cov = dist.cov,
                     trap.cov = data.frame(trap = 1:nrow(ihd$traps), forest_volumn = "x", brand = "a")),
    "Covariates 'noise', 'forest_volumn' in 'loc.cov' are ignored in favour of the covariates attached"),
    "Covariate 'noise' in 'dist.cov'"), "Covariate 'forest_volumn' in 'trap.cov'")
  cov <- dat$par.extend$data$mask
  expect_equal(cov$noise, mask$noise)
  expect_equal(cov$forest_volumn, mask$forest_volumn)
  expect_true(all(c("elevation", "road") %in% colnames(cov)))
  expect_equal(colnames(dat$par.extend$data$trap), c("trap", "brand"))
  expect_equal(names(dat$arg.input$dist.cov), "road")
  expect_warning(read.acre(captures = ihd$capt, traps = ihd$traps, mask = mask,
                           session.cov = data.frame(session = 1, noise = 1)),
                 "Covariate 'noise' in 'session.cov'")
})

test_that("covariates attached to masks work for several sessions", {
  d <- mul_ses
  masks <- create.mask(d$traps, buffer = d$control.mask$buffer, nx = 10)
  masks_cov <- lapply(seq_along(masks), function(s){
    m <- data.frame(masks[[s]], habitat = s*10 + seq_len(nrow(masks[[s]])))
    attr(m, "area") <- attr(masks[[s]], "area")
    m
  })
  dat <- read.acre(captures = d$capt, traps = d$traps, mask = masks_cov)
  cov <- dat$par.extend$data$mask
  for(s in 1:2){
    expect_equal(cov$habitat[cov$session == s], masks_cov[[s]]$habitat)
    expect_equal(cov$mask[cov$session == s], seq_len(nrow(masks[[s]])))
  }
  masks_bad <- masks_cov
  names(masks_bad[[2]])[3] <- "other"
  expect_error(read.acre(captures = d$capt, traps = d$traps, mask = masks_bad), "same covariates attached")
})

test_that("fit.acre() uses covariates attached to a new mask, and needs them", {
  mask <- ihd.mask.cov()
  dat <- read.acre(captures = ihd$capt, traps = ihd$traps, mask = mask)
  model <- list(D = ~ noise)
  fix <- list(fix = list(g0 = 1))
  ## A new mask without the attached covariates.
  expect_error(fit.acre(dat, model = model, control.mask = list(buffer = 30, nx = 10), tracing = FALSE),
               "'noise', 'forest_volumn' were attached to the mask given to read.acre\\(\\)")
  ## A new mask with them attached gives the same fit as providing it to read.acre().
  mask_new <- mask[seq(1, nrow(mask), by = 2), ]
  attr(mask_new, "area") <- 2*attr(mask, "area")
  fit_new <- fit.quietly(dat, list(model = model, mask = mask_new, optim.opts = fix))
  fit_read <- fit.quietly(read.acre(captures = ihd$capt, traps = ihd$traps, mask = mask_new),
                          list(model = model, optim.opts = fix))
  expect_equal(get.par.extend.data(fit_new)$mask$noise, mask_new$noise)
  expect_equal(fit_new$loglik, fit_read$loglik)
  expect_equal(coef(fit_new), coef(fit_read))
  ## Attached to a new mask, in favour of interpolated covariates.
  dat_loc <- read.acre(captures = ihd$capt, traps = ihd$traps, loc.cov = ihd$loc.cov,
                       control.mask = list(buffer = ihd$control.mask$buffer, nx = 16))
  expect_warning(utils::capture.output(
    fit_att <- fit.acre(dat_loc, model = model, mask = mask_new[, c("x", "y", "noise")], optim.opts = fix,
                        tracing = FALSE)),
    "Covariate 'noise' in 'loc.cov' is ignored")
  expect_equal(get.par.extend.data(fit_att)$mask$noise, mask_new$noise)
})

test_that("covariates attached to masks are checked", {
  mask <- ihd.mask.cov()
  bad <- mask
  names(bad)[1:2] <- c("east", "north")
  expect_error(read.acre(captures = ihd$capt, traps = ihd$traps, mask = bad), "columns 'x' and 'y'")
  bad <- mask
  bad$session <- 1
  expect_error(read.acre(captures = ihd$capt, traps = ihd$traps, mask = bad), "cannot be called 'session' or 'mask'")
  ## Missing values are not interpolated, in read.acre() or fit.acre().
  bad <- mask
  bad$noise[c(3, 7)] <- NA
  bad$forest_volumn[5] <- NA
  msg <- "cannot have missing values, but 'noise' has 2, 'forest_volumn' has 1 for session 1. Use 'loc.cov'"
  expect_error(read.acre(captures = ihd$capt, traps = ihd$traps, mask = bad), msg, fixed = TRUE)
  dat <- read.acre(captures = ihd$capt, traps = ihd$traps, mask = mask)
  expect_error(fit.acre(dat, model = list(D = ~ noise), mask = bad, tracing = FALSE), msg, fixed = TRUE)
})

test_that("create.mask() returns data frames that covariates can be attached to", {
  mask <- create.mask(ihd$traps, buffer = ihd$control.mask$buffer, nx = 16)
  expect_true(is.data.frame(mask))
  expect_equal(colnames(mask), c("x", "y"))
  expect_equal(attr(mask, "buffer"), ihd$control.mask$buffer)
  expect_true(attr(mask, "area") > 0)
  masks <- create.mask(mul_ses$traps, buffer = mul_ses$control.mask$buffer, nx = 10)
  expect_true(is.list(masks) && all(vapply(masks, is.data.frame, logical(1))))
  ## The same as the mask read.acre() creates, which is stored as a matrix.
  dat <- read.acre(captures = ihd$capt, traps = ihd$traps, control.mask = list(buffer = ihd$control.mask$buffer, nx = 16))
  expect_true(is.matrix(dat$mask[[1]]))
  expect_equal(dat$mask[[1]], as.matrix(mask), ignore_attr = TRUE)
  expect_equal(attr(dat$mask[[1]], "area"), attr(mask, "area"))
  ## Masks provided as data frames are stored as matrices, with or without covariates attached.
  expect_equal(read.acre(captures = ihd$capt, traps = ihd$traps, mask = mask)$mask, dat$mask)
  mask$noise <- mask$x
  dat_cov <- read.acre(captures = ihd$capt, traps = ihd$traps, mask = mask)
  expect_equal(dat_cov$mask, dat$mask)
  fit <- fit.quietly(dat, list(model = list(D = ~ noise), mask = mask, optim.opts = list(fix = list(g0 = 1))))
  expect_true(is.matrix(get.mask(fit)[[1]]))
  expect_equal(get.mask(fit), dat$mask)
  ## A mask with only two unnamed columns is still allowed.
  m <- unname(as.matrix(mask[, 1:2]))
  attr(m, "area") <- attr(mask, "area")
  expect_equal(colnames(read.acre(captures = ihd$capt, traps = ihd$traps, mask = m)$mask[[1]]), c("x", "y"))
})
