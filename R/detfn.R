## The argument dx and each component of det.par should be vectors of the same length;
## of course, if any of them is a scalar, there is no problem.

## The argument orientation is used for ss.dir, which is not supported yet.

det.prob <- function(det.fn, det.par, dx, ss.link = NULL, orientation = NULL) {
  g0 <- NULL
  sigma <- NULL
  lambda0 <- NULL
  z <- NULL
  shape <- NULL
  shape.1 <- NULL
  shape.2 <- NULL
  b0.ss <- NULL
  b1.ss <- NULL
  
  
  for(i in names(det.par)) assign(i, det.par[[i]])
  
  if(det.fn == 'hn'){
    out <- g0 * exp((-1 * dx^2) / (2 * sigma^2))
  } else if(det.fn == 'hhn'){
    out <- 1 - exp((-1) * lambda0 * exp(dx^2 * (-0.5) / sigma^2))
  } else if(det.fn == 'hr'){
    out <- g0 * (1 - exp((-1) * (dx/sigma)^(-z)))
  } else if(det.fn == 'th'){
    out <- 0.5 - 0.5 * erf(dx / scale - shape)
  } else if(det.fn == 'lth'){
    out <- 0.5 - 0.5 * erf(shape.1 - exp(shape.2 - scale * dx))
  } else if(det.fn == 'ss'){
    ## The ss model is a little different: if the simulated ss > cutoff, the detection probability is 1,
    ## and 0 otherwise, so we return E(SS) instead of the detection probability.
    if(ss.link == 'identity'){
      out <- b0.ss - b1.ss * dx
    } else if(ss.link == 'log'){
      out <- exp(b0.ss - b1.ss * dx)
    } else if(ss.link == 'spherical'){
      out <- ifelse(dx > 1, b0.ss - 20 * log10(dx) - b1.ss * (dx - 1), b0.ss)
    }
  }
  
  return(out)
  
}



