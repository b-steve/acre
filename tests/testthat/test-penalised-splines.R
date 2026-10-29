## Tests for smooth terms in the density model, with and without penalisation.

## The ihd demo data on a coarse mask.
ihd.coarse <- function(nx = 16){
  read.acre(captures = ihd$capt, traps = ihd$traps,
            control.mask = list(buffer = ihd$control.mask$buffer, nx = nx), loc.cov = ihd$loc.cov)
}

test_that("smooth design matrices for the fitting data match mgcv's", {
  set.seed(1)
  d <- data.frame(x = runif(200), y = runif(200), z = runif(200), f = factor(sample(c("a", "b"), 200, TRUE)))
  d$gam.resp <- 1
  G <- mgcv::gam(gam.resp ~ f + s(x, y, k = 12) + s(z, k = 5, fx = TRUE), data = d, fit = FALSE)
  ## Unpenalised, the columns are mgcv's design matrix, in the same order.
  info <- smooth.info(G, d, penalise = FALSE)
  des <- smooth.design(info, d)
  expect_equal(unname(des$X), unname(G$X))
  expect_equal(ncol(des$Z), 0)
  ## Penalised, the columns are those from smooth2random(); the fx = TRUE smooth is not penalised.
  info <- smooth.info(G, d, penalise = TRUE)
  des <- smooth.design(info, d)
  sm <- G$smooth[[1]]
  sm$X <- G$X[, sm$first.para:sm$last.para]
  re <- mgcv::smooth2random(sm, "", type = 2)
  expect_equal(unname(des$Z), unname(re$rand[[1]]), ignore_attr = TRUE)
  expect_equal(unname(des$X[, 3:4]), unname(re$Xf))
  expect_equal(unname(des$X[, -(3:4)]), unname(G$X[, c(1:2, 14:17)]))
  expect_equal(des$re_fac, rep(1, ncol(des$Z)))
  expect_equal(des$sd_labels, "s(x,y)")
  ## Building the matrices for new data does not set the model up again.
  new <- d[c(5, 1, 7), ]
  expect_equal(smooth.design(info, new)$Z, des$Z[c(5, 1, 7), ], ignore_attr = TRUE)
  expect_equal(smooth.design(info, new)$X, des$X[c(5, 1, 7), ], ignore_attr = TRUE)
})

test_that("predicted density surfaces match the fitted density for smooth models", {
  dat <- ihd.coarse()
  ## Unpenalised and penalised smooths.
  for(model in list(~ s(x, y, k = 8, fx = TRUE), ~ s(x, y, k = 8))){
    fit <- fit.quietly(dat, list(model = list(D = model), detfn = 'hn',
                                 optim.opts = list(fix = list(g0 = 1))))
    pred <- D.predict.for.plot(fit, se.fit = TRUE)
    expect_equal(pred$est, unname(fit$D.mask[[1]]), tolerance = 1e-8)
    expect_true(all(is.finite(pred$std) & pred$std > 0))
  }
})

test_that("the penalised likelihood is the Laplace approximation to the marginal likelihood", {
  dat <- ihd.coarse()
  fit <- fit.quietly(dat, list(model = list(D = ~ y + s(x, k = 6)), detfn = 'hn',
                               optim.opts = list(fix = list(g0 = 1))))
  expect_s3_class(fit, 'acre')
  expect_equal(fit$D.smooth$term, "s(x)")
  ## The fitted coefficients for D, on the scale used for fitting, are the fixed coefficients,
  ## then the random effects.
  est <- fit$output.tmb$tmb_output_og$est
  n_u <- sum(fit$D.smooth$n.re)
  D_all <- est[names(est) == 'D']
  n_b <- length(D_all) - n_u
  beta <- D_all[1:n_b]
  X <- cbind(1, fit$output.tmb$DX_standardized$DX_mask$D)
  sigma <- exp(est[names(est) == 'sigma'])
  A <- fit$output.tmb$area_unit
  ## The likelihood for the call model with a halfnormal detection function.
  mask <- get.mask(fit)[[1]]
  traps <- as.matrix(ihd$traps)
  dists <- as.matrix(dist(rbind(mask, traps)))[1:nrow(mask), nrow(mask) + 1:nrow(traps)]
  p <- exp(-dists^2 / (2 * sigma^2))
  p_dot <- 1 - apply(1 - p, 1, prod)
  capt <- fit$output.tmb$data.full
  capt <- capt[!is.na(capt$ID), ]
  W <- matrix(0, length(unique(capt$ID)), nrow(traps))
  W[cbind(match(capt$ID, unique(capt$ID)), capt$trap)] <- capt$bincapt
  f_w <- exp(W %*% t(log(p)) + (1 - W) %*% t(log(1 - p)))
  joint <- function(u){
    D <- exp(drop(X %*% c(beta, u)))
    -A * sum(D * p_dot) + sum(log(A * f_w %*% D)) + sum(dnorm(u, 0, fit$D.smooth$sd, log = TRUE))
  }
  o <- stats::optim(D_all[-(1:n_b)], function(u) -joint(u), method = 'BFGS', control = list(reltol = 1e-14))
  H <- stats::optimHess(o$par, function(u) -joint(u))
  laplace <- -o$value + n_u / 2 * log(2 * pi) - 0.5 * as.numeric(determinant(H)$modulus)
  expect_equal(fit$loglik, laplace, tolerance = 1e-6)
  expect_equal(unname(o$par), unname(D_all[-(1:n_b)]), tolerance = 1e-5)
})

test_that("penalised fits report random effects separately from the parameters", {
  dat <- ihd.coarse()
  fit <- fit.quietly(dat, list(model = list(D = ~ s(x, y, k = 8)), detfn = 'hn',
                               optim.opts = list(fix = list(g0 = 1))))
  n_u <- sum(fit$D.smooth$n.re)
  re_pattern <- "\\.r[0-9]+(_link)?$"
  ## The random effects are not among the parameters.
  for(type in c('linked', 'fitted')){
    expect_false(any(grepl(re_pattern, names(coef(fit, types = type)))))
    expect_false(any(grepl(re_pattern, names(stdEr(fit, types = type)))))
    expect_false(any(grepl(re_pattern, rownames(confint(fit, types = type)))))
    expect_false(any(grepl(re_pattern, rownames(vcov(fit, types = type)))))
  }
  n_D <- length(coef(fit, types = 'linked', pars = 'D'))
  expect_equal(dim(vcov(fit, types = 'linked', pars = 'D')), c(n_D, n_D))
  ## They are reported separately instead, and used internally with D's coefficients.
  expect_equal(nrow(fit$D.re), n_u)
  expect_equal(fit$D.re$term, rep("s(x,y)", n_u))
  expect_true(all(fit$D.re$se > 0))
  all_D <- coef(fit, types = 'linked', pars = 'D', re = TRUE)
  expect_equal(unname(all_D[-(1:n_D)]), fit$D.re$est)
  ## The random effects are not parameters for AIC, but their standard deviation is.
  expect_equal(AIC(fit)$npar, length(coef(fit)) + 1)
  s <- summary(fit)
  expect_false(any(grepl(re_pattern, names(s$coefs))))
  expect_output(print(s), "Penalised smooths for D")
  ## Bootstrap results record the random effects, so bootstrap density surfaces include them.
  b <- suppressWarnings(boot.acre(fit, N = 3, seed = 1))
  expect_equal(sum(grepl(re_pattern, colnames(b$boot$boots))), n_u)
  expect_false(any(grepl(re_pattern, names(coef(b, types = 'linked', correct.bias = TRUE)))))
  expect_false(any(grepl(re_pattern, rownames(confint(b)))))
  expect_false(any(grepl(re_pattern, names(stdEr(b)))))
  p <- D.predict.for.plot(b, se.fit = TRUE)
  expect_equal(p$est, unname(fit$D.mask[[1]]), tolerance = 1e-8)
})

test_that("penalised smooths are checked", {
  dat <- ihd.coarse(nx = 8)
  expect_error(fit.acre(dat, model = list(D = ~ s(x, k = 5), sigma = ~ s(x, k = 5)), detfn = 'hn', tracing = FALSE),
               "not available for parameters other than D")
  expect_error(fit.acre(dat, model = list(D = ~ s(x, k = 5)), detfn = 'hn', tracing = FALSE,
                        optim.opts = list(gr.skip = TRUE)),
               "'gr.skip' cannot be TRUE")
})
