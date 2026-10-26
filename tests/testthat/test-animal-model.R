## Likelihood properties of the animal model (Stevenson et al., 2021) that should hold exactly,
## so these tests compare fits with each other rather than with stored values.

## One session of the ind_bearing_dist data, with homogeneous density, on a coarse mask. The buffer
## is wider than the demo's, so that detection probability falls off within the mask and the ESA
## depends on the detection parameters; otherwise the conditional likelihood test below is insensitive.
## One detection has a missing bearing, which should simply contribute no bearing information.
animal.data <- function(rotate = FALSE){
  d <- ind_bearing_dist
  capt <- d$capt[d$capt$session == 1, c("session", "ID", "trap", "animal_ID", "bearing", "dist")]
  capt$bearing[3] <- NA
  traps <- as.matrix(d$traps)
  mask <- create.mask(traps, buffer = 40, nx = 16)
  if(rotate){
    ## Rotate everything 90 degrees clockwise: (x, y) -> (y, -x). Bearings are clockwise from
    ## north, so they all increase by pi/2.
    rot <- function(m) cbind(m[, 2], -m[, 1])
    traps <- rot(traps)
    mask <- structure(rot(mask), area = attr(mask, "area"), buffer = attr(mask, "buffer"))
    capt$bearing <- (capt$bearing + pi/2) %% (2*pi)
  }
  colnames(traps) <- c("x", "y")
  colnames(mask) <- c("x", "y")
  suppressWarnings(read.acre(capt, traps, mask = mask, survey.length = d$survey.length[1]))
}

fit.animal <- fit.quietly(animal.data(), list(detfn = "hn"))

test_that("animal model fits are invariant to rotating the survey", {
  ## Undetected traps carry placeholder bearings of zero, which do not rotate, and the missing
  ## bearing would become zero if it were not kept as NA; either would change the rotated fit.
  fit.rot <- fit.quietly(animal.data(rotate = TRUE), list(detfn = "hn"))
  expect_equal(coef(fit.rot, types = 'linked'), coef(fit.animal, types = 'linked'), tolerance = 1e-5)
})

test_that("conditional likelihood matches the full likelihood for one session", {
  ## With one session and homogeneous density, the full likelihood factorises as
  ## Pois(n; D * ESA) * f(captures | n), so both give the same estimates of the other parameters,
  ## and n / ESA from the conditional fit equals the full-likelihood estimate of D.
  fit.cl <- fit.quietly(animal.data(), list(detfn = "hn", CL = TRUE))
  pars <- setdiff(names(coef(fit.animal, types = 'linked')), "D_link")
  expect_equal(coef(fit.cl, types = 'linked')[pars], coef(fit.animal, types = 'linked')[pars], tolerance = 1e-4)
  n.animals <- length(unique(ind_bearing_dist$capt$animal_ID[ind_bearing_dist$capt$session == 1]))
  expect_equal(unname(n.animals / esa(fit.cl)$value), unname(coef(fit.animal, types = 'fitted')['D']), tolerance = 1e-4)
})

test_that("conditional likelihood matches the full likelihood for the call model", {
  d <- simple_hhn
  dat <- read.acre(d$capt, d$traps, control.mask = list(buffer = 30, nx = 15))
  fit <- fit.quietly(dat, list(detfn = d$detfn))
  fit.cl <- fit.quietly(dat, list(detfn = d$detfn, CL = TRUE))
  pars <- setdiff(names(coef(fit, types = 'linked')), "D_link")
  expect_equal(coef(fit.cl, types = 'linked')[pars], coef(fit, types = 'linked')[pars], tolerance = 1e-4)
  n.calls <- length(unique(d$capt$ID))
  expect_equal(unname(n.calls / esa(fit.cl)$value), unname(coef(fit, types = 'fitted')['D']), tolerance = 1e-4)
})
