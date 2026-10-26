#the comprehensive tests on the default (fine) masks take several minutes, so they live in tests/full
#and are run with testthat::test_dir("tests/full", load_package = "source")

#fit a demo model the same way as demo_fit(), but on a coarse mask (nx points across), which is
#much faster; the estimates differ from the ones on the default mask, but are just as reproducible
fit_demo_coarse = function(data_name, nx = 20, gradient_free = FALSE){
  dat = get(data_name)
  if(!is.null(dat$control_create_mask)) dat$control.mask = dat$control_create_mask
  dat$control.mask$nx = nx
  fit_args = list(model = dat$model, detfn = dat$detfn, ss.opts = dat$ss.opts,
                  optim.opts = list(sv = dat$sv, fix = dat$fix, gr.skip = gradient_free))
  dat[c('model', 'sv', 'fix', 'bounds', 'detfn', 'ss.opts', 'control.capt', 'control_create_mask')] = NULL
  fit_quietly(do.call('read.acre', dat), fit_args)
}

fit_quietly = function(dat, fit_args = list()){
  out = NULL
  utils::capture.output(out <- suppressMessages(suppressWarnings(
    do.call('fit.acre', c(list(data = dat, tracing = FALSE), fit_args)))))
  out
}
