# Buffered trace printer for TMB objectives (paired fn/grad)

Returns a small tap object with callbacks that you attach to the
objective and gradient produced by
[`TMB::MakeADFun()`](https://rdrr.io/pkg/TMB/man/MakeADFun.html). While
optimizing (e.g. via
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html)), it prints a
right-aligned, fixed-width row **only when** the objective value and
gradient were both evaluated at the **same** parameter vector. This
guarantees the `mgc` column (maximum absolute gradient component)
corresponds to the printed `fval`.

## Usage

``` r
make_buffer_printer(
  trace_cols,
  only_improvements = TRUE,
  step = c("euclid", "max", "none"),
  show_mgc = TRUE,
  mgc_tol = 1e-10
)
```

## Arguments

- trace_cols:

  Character vector of parameter names to print as columns (typically
  `names(obj$par)` or a subset).

- only_improvements:

  Logical; if `TRUE` (default) print only when `fval` strictly improves
  over the best so far (an “outer-like” trace). Set `FALSE` to print
  every paired evaluation.

- step:

  One of `"euclid"`, `"max"`, `"none"`. Controls the optional step-size
  column between successive *printed* parameter vectors.

- show_mgc:

  Logical; if `TRUE`, include a `mgc` column equal to `max(abs(grad))`
  at the matched point. Set `FALSE` if you do not have a gradient (e.g.
  `gr.skip = TRUE`).

- mgc_tol:

  Numeric tolerance used to decide whether `fn` and `gr` parameter
  vectors “match” (default `1e-10`).

## Value

A list with callbacks:

- `on_fn(par_vec, fval, par_names, G_matrix_list)`: record an objective
  evaluation.

- `on_gr(par_vec, grad_vec, par_names, G_matrix_list)`: record a
  gradient evaluation.

Each callback returns `invisible(NULL)`; printing happens as a side
effect.

## Details

Because optimizers often evaluate the objective and gradient at
different points and in different orders (e.g. line search), the tap
uses an internal buffer to: (a) remember the most recent `fn` and `gr`
evaluations, and (b) print a row only when those evaluations match in
parameters (within `mgc_tol`).

Columns are `iter` (count of printed rows), `fval`, optional `mgc`,
optional `step` (Euclidean or max coordinate change since the previous
*printed* row), followed by the parameter columns in `trace_cols`.

## Usage


    tap <- make_buffer_printer(trace_cols = names(obj$par),
                               step = "euclid",
                               show_mgc = TRUE)
    obj$fn <- with_fn_tap(obj$fn, tap, par_names = names(obj$par))
    obj$gr <- with_gr_tap(obj$gr, tap, par_names = names(obj$par))

    opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                                   control = list(trace = 0))

## See also

[`TMB::MakeADFun()`](https://rdrr.io/pkg/TMB/man/MakeADFun.html),
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html)
