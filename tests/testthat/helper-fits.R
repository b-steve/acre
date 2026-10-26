## Fit a demo model the same way as demo.fit(), but on a coarse mask (nx points across), which is
## much faster; the estimates differ from the ones on the default mask, but are just as reproducible.
fit.demo.coarse <- function(data.name, nx = 16, gradient.free = FALSE){
  dat <- get(data.name)
  if(!is.null(dat$control_create_mask)) dat$control.mask <- dat$control_create_mask
  dat$control.mask$nx <- nx
  fit_args <- list(model = dat$model, detfn = dat$detfn, ss.opts = dat$ss.opts,
                  optim.opts = list(sv = dat$sv, fix = dat$fix, gr.skip = gradient.free))
  dat[c('model', 'sv', 'fix', 'bounds', 'detfn', 'ss.opts', 'control.capt', 'control_create_mask')] <- NULL
  fit.quietly(do.call('read.acre', dat), fit_args)
}

fit.quietly <- function(dat, fit.args = list()){
  out <- NULL
  utils::capture.output(out <- suppressMessages(suppressWarnings(
    do.call('fit.acre', c(list(data = dat, tracing = FALSE), fit.args)))))
  out
}
